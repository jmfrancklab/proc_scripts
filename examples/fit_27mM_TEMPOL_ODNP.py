r"""High-concentration TEMPOL ODNP: T₁ map, k\ :sub:`ρ`, and motional narrowing.

Goal
----
Map T₁(p) for high-concentration hydroxytempo in water, then fit both
R\ :sub:`1,w`\(p) and k\ :sub:`ρ`\(p) to an Arrhenius motional-narrowing
model.

Physics recap
-------------
From Franck & Han (2019) / DCCT paper:

    M₀E(p) = M₀ − M₀ A s(p) / R₁(p)

where A = k\ :sub:`σ` s\ :sub:`max` C\ :sub:`SL` |ω\ :sub:`e`/ω\ :sub:`H`|
[units s⁻¹], s(p) = p/(p+p\ :sub:`½`), and

    R₁(p) = R₁,w(p) + k\ :sub:`ρ` C\ :sub:`SL`
           = [T₁,w,0 + ΔT₁,w p]⁻¹ + C / (a₀ + a₁ p)

with [a₀, a₁] the linear fit to k\ :sub:`ρ`\ ⁻¹(p).

Motional narrowing (extreme-narrowing limit, ω τ\ :sub:`c` ≪ 1 at 15 MHz):

    quantity(p) = offset + scale · τ\ :sub:`c`\(p)
    τ\ :sub:`c`\(p) = τ₀ exp[E\ :sub:`a`/R · (1/T(p) − 1/T₀)]
    T(p) = T₀ + α p

The motional narrowing fit is **only attempted when there are ≥ 4 FIR power
points**. With the example two-point dataset the script produces valid T₁
maps and the k\ :sub:`ρ`\ ⁻¹ linear fit, but skips the Arrhenius fit.

Key proc_scripts conventions followed
--------------------------------------
* ``acq_params`` is accessed via ``s.get_prop("acq_params")[key]`` after
  ``find_file``, NOT via ``h5py`` group attributes directly.
* ``FIR_rep`` (µs, no ``_us`` suffix) is the correct key for the FIR
  repetition time, matching ``proc_FIR.py`` in the proc_scripts examples.
* The axis swap after ``rough_table_of_integrals`` follows the exact pattern
  in ``proc_Ep.py``: restore the full structured ``orig_axis`` first, then
  extract the ``"power"`` field — never set ``set_error("indirect", ...)``
  after the axis has already been renamed to ``"power"``.
* Normalisation (``s /= s["indirect", 0:1]``) and a synthetic error
  (``s.set_error(0.01 * s["indirect", 0].item())``) are applied before the
  axis swap, exactly as in ``proc_Ep.py``.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import re
import sympy as sp
import h5py
from pyspecProcScripts.generate_coordinates_from_log import (
    generate_coordinates_from_log,
)

# {{{ user-editable block
thisfile, thisexptype, nodename = (
    "260526_hydroxytempo_ODNP_2.h5",
    "B27/ODNP",
    "ODNP",
)

# Motional-narrowing Arrhenius parameters (fixed, not fit).
# Change if the solvent or radical differ substantially from TEMPOL/water.
T0 = 298.15  # K   — ambient/reference temperature
heating_K_per_W = 5.0  # K/W — microwave-induced sample heating
Ea = 15e3  # J/mol — activation energy for τ_c
tau0 = 30e-12  # s   — pre-exponential factor for τ_c
MIN_FIR_FOR_MN = 4  # min FIR power points needed for motional narrowing fit
# }}}


# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────
def node_power(nodename):
    """Return microwave power in watts from a FIR node name."""
    m = re.search(r"([0-9]+(?:p[0-9]+)?)dBm", nodename)
    if m:
        return prscr.dBm2power(float(m.group(1).replace("p", "."))).item()
    if "noPower" in nodename:
        return 0.0
    raise ValueError(f"Cannot infer power from node name {nodename!r}")


def sorted_fir_nodes(h5file, exp_type):
    """Return list of (nodename, power_W) for all FIR_ nodes, ascending."""
    filename = psd.search_filename(h5file, exp_type=exp_type, unique=True)
    with h5py.File(filename, "r") as f:
        nodes = [k for k in f.keys() if k.startswith("FIR_")]
    return sorted([(n, node_power(n)) for n in nodes], key=lambda x: x[1])


def fir_signal_range(s, direct="t2", peak_lower_thresh=0.5):
    signal_pathway = s.get_prop("coherence_pathway")
    pathway_data = prscr.select_pathway(s.C, signal_pathway)
    frq_center, frq_half = prscr.find_peakrange(
        pathway_data, direct=direct, peak_lower_thresh=peak_lower_thresh
    )
    return (frq_center - frq_half, frq_center + frq_half)


# ─────────────────────────────────────────────────────────────────────────────
# pre-flight
# ─────────────────────────────────────────────────────────────────────────────
FIR_nodes = sorted_fir_nodes(thisfile, thisexptype)
n_fir = len(FIR_nodes)
print(f"Found {n_fir} FIR node(s): {[n for n, _ in FIR_nodes]}")
if n_fir < MIN_FIR_FOR_MN:
    print(
        f"WARNING: only {n_fir} FIR power point(s). "
        f"Motional narrowing fit requires ≥ {MIN_FIR_FOR_MN} points — "
        f"will be skipped. k_ρ⁻¹ linear fit is exact (zero residual) "
        f"through 2 points; treat as interpolation only."
    )

# ─────────────────────────────────────────────────────────────────────────────
# single figlist_var wraps all figures — post_proc pattern
# ─────────────────────────────────────────────────────────────────────────────
with psd.figlist_var() as fl:
    fl.basename = thisfile

    # {{{ load and integrate the progressive-enhancement (Ep) dataset
    #     following proc_Ep.py exactly
    s = psd.find_file(
        thisfile,
        exp_type=thisexptype,
        expno=nodename,
        lookup=prscr.lookup_table,
        fl=fl,
    )
    if s.get_prop("log") is None:
        s = prscr.attach_log_data_from_file(s, thisfile, thisexptype)
    m = re.search(r".*ODNP.*v([0-9]+)$", s.get_prop("postproc_type"))
    if m:
        vernum = int(m.groups()[0])
    else:
        raise IOError(
            f"Unexpected postproc_type: {s.get_prop('postproc_type')!r}"
        )
    if vernum < 6:
        s = generate_coordinates_from_log(s, fl=fl)
    # stash structured indirect axis before rough_table collapses it
    orig_axis = s["indirect"]
    orig_axis_error = s.get_error("indirect")
    s["indirect"] = s["indirect"]["time"]
    s.set_units("indirect", "s")
    s.set_error("indirect", orig_axis_error["time"])
    s, _ = prscr.rough_table_of_integrals(s, fl=fl)
    # synthetic error — proc_Ep.py: 1% of first point so relative error bars visible
    s.set_error(s["indirect", 0].item() * 0.01)
    # normalise to thermal (p=0) integral
    s /= s["indirect", 0:1]
    # restore full structured axis, then extract power field
    # (this is the exact proc_Ep.py pattern; set_error BEFORE rename)
    s["indirect"] = orig_axis
    s.set_error("indirect", orig_axis_error)
    s.set_error("indirect", s.get_error("indirect")["power"])
    s["indirect"] = s["indirect"]["power"]
    s.set_units("indirect", "W").rename("indirect", "power")
    Ep = s
    # get acq_params via get_prop (pyspecdata standard, not h5py direct)
    acq = Ep.get_prop("acq_params")
    C = float(acq["concentration"])
    chemical = acq.get("chemical", b"TEMPOL")
    if isinstance(chemical, bytes):
        chemical = chemical.decode()
    sample_label = f"{1e3 * C:g} mM {chemical}"
    # }}}

    # {{{ fit inversion-recovery (FIR) experiments → R₁(p)
    #     following proc_FIR.py, with FIR repetition key fallback for older files
    Mi, R1, vd = sp.symbols("M_inf R_1 vd", real=True)
    R1p = psd.ndshape([("power", n_fir)]).alloc(dtype=np.float64)
    R1p.setaxis("power", [pw for _, pw in FIR_nodes]).set_units("power", "W")

    for j, (thisnodename, _) in enumerate(FIR_nodes):
        fl.basename = thisnodename
        s = psd.find_file(
            thisfile,
            exp_type=thisexptype,
            expno=thisnodename,
            lookup=prscr.lookup_table,
        )
        if "nScans" in s.dimlabels:
            s = prscr.clock_correct(s)
        s = s.squeeze()
        signal_range = fir_signal_range(s)
        s, ax_last = prscr.rough_table_of_integrals(
            s, fl=fl, signal_range=signal_range
        )
        this_acq = s.get_prop("acq_params")
        # FIR repetition key changed across file versions
        FIR_rep = this_acq.get("FIR_rep", this_acq["FIR_rep_us"])
        W = FIR_rep * 1e-6 + this_acq["acq_time_ms"] * 1e-3
        # scale R₁ guess to vd axis units (guard for missing units)
        vd_units = s.get_units("vd")
        prefactor_scaling = (
            10 ** psd.det_unit_prefactor(vd_units)
            if vd_units is not None
            else 1.0
        )
        s = psd.lmfitdata(s)
        s.functional_form = Mi * (1 - (2 - sp.exp(-W * R1)) * sp.exp(-vd * R1))
        s.set_guess(
            M_inf=dict(
                value=s.max().item(),
                min=0.1 * s.max().item(),
                max=1.5 * s.max().item(),
            ),
            R_1=dict(
                value=0.8 * prefactor_scaling,
                min=0.01 * prefactor_scaling,
                max=100 * prefactor_scaling,
            ),
        )
        s.fit()
        R1p["power", j] = s.output("R_1")
        s_fit = s.eval(200)
        psd.plot(s_fit, ax=ax_last, alpha=0.5)
        ax_last.text(
            0.5,
            0.5,
            f"{thisnodename} RESULT: %s" % s.latex(),
            ha="center",
            va="center",
            color=s_fit.get_plot_color(),
            transform=ax_last.transAxes,
        )
    fl.basename = thisfile  # restore after FIR loop
    # }}}

    # {{{ flip index (ascending → descending power in Ep)
    flip_idx = np.argmax(Ep.getaxis("power")) + 1
    # }}}

    # {{{ T₁,w(p) polynomial and k_ρ⁻¹(p) linear fit
    T10_p = np.r_[
        acq["T1water_cold"],
        (acq["T1water_hot"] - acq["T1water_cold"]) / acq["max_power"],
    ]
    R10_p = 1.0 / R1p.fromaxis("power").eval_poly(T10_p, "power")
    krho = (R1p - R10_p) / C
    krho_inv = C / (R1p - R10_p)
    krho_inv_coeff = krho_inv.polyfit("power", order=1)
    # }}}

    # {{{ symbolic R₁(p) expression with frozen k_ρ⁻¹ linear coefficients
    M0, A, phalf, p = sp.symbols("M0 A phalf power", real=True)
    R1p_expr = (T10_p[0] + T10_p[1] * p) ** -1 + (
        C / (krho_inv_coeff[0] + krho_inv_coeff[1] * p)
    )
    p_max = max(Ep.getaxis("power").max(), R1p.getaxis("power").max())
    powers_fine = psd.nddata(np.r_[0:p_max:300j], "power").set_units(
        "power", "W"
    )
    R1p_fit = powers_fine.fromaxis("power").run(
        sp.lambdify(p, R1p_expr, "numpy")
    )
    krho_fit = powers_fine.fromaxis("power").run(
        lambda x: 1.0 / (krho_inv_coeff[0] + krho_inv_coeff[1] * x)
    )
    # }}}

    # {{{ motional narrowing fits (only if ≥ MIN_FIR_FOR_MN points)
    do_mn_fit = n_fir >= MIN_FIR_FOR_MN
    if do_mn_fit:
        offset_sym, scale_sym = sp.symbols("offset scale", real=True)
        mn_tau_expr = tau0 * sp.exp(
            Ea / 8.314 * (1.0 / (T0 + heating_K_per_W * p) - 1.0 / T0)
        )
        mn_form = offset_sym + scale_sym * mn_tau_expr

        def fit_motional_narrowing(y_nddata, p_fine, guess_scale):
            f = psd.lmfitdata(y_nddata.real)
            f.functional_form = mn_form
            f.set_guess(
                offset=dict(value=float(y_nddata.real.data.min()), vary=True),
                scale=dict(value=guess_scale, vary=True),
            )
            f.fit()
            return f.eval(p_fine)

        R10p_for_mn = 1.0 / R1p.fromaxis("power").eval_poly(T10_p, "power")
        R10p_mn_fit = fit_motional_narrowing(R10p_for_mn, powers_fine, 1e10)
        krho_mn_fit = fit_motional_narrowing(krho, powers_fine, 1e13)
    # }}}

    # {{{ fit M₀E(p) (normalised integrals) vs microwave power
    # After normalisation Ep ≈ E(p) = 1 − ε(p), so M0 ≈ 1.
    # Functional form: E(p) = M0 − M0 A s(p) / R1(p)
    # A = kσ smax CSL |ωe/ωH|  [s⁻¹], s(p) = p/(p+phalf)
    sp_expr = p / (p + phalf)
    Ep_fit = psd.lmfitdata(Ep["power", :flip_idx].real)
    Ep_fit.functional_form = M0 - (M0 * A * sp_expr) / R1p_expr
    # A guess: A ≈ (1 − E(p_flip)) · R₁(0)   (s ≈ 1 at high power)
    A_guess = 1.0 - (
        R1p["power", 0].real.item()
        * (Ep["power", flip_idx].real.item() / Ep["power", 0].real.item())
    )
    M0_val = Ep["power", 0].real.item()
    Ep_fit.set_guess(
        M0=dict(value=M0_val, min=0.1 * M0_val, max=2.0 * M0_val),
        A=dict(value=A_guess, min=0.2 * A_guess, max=3.0 * A_guess),
        phalf=dict(
            value=float(acq.get("guessed_phalf", 0.2)), min=0.05, max=1.0
        ),
    )
    Ep_fit.fit()
    Ep_fit_curve = Ep_fit.eval(100)
    ksig = Ep_fit.output("A") * float(acq["guessed_MHz_to_GHz"]) * 1e-3 / C
    # }}}

    # ─── plots ────────────────────────────────────────────────────────────────

    # {{{ normalised E(p) enhancement
    fl.next(f"{sample_label} – E(p)")
    Ep_progressive = Ep["power", :flip_idx].C.set_plot_color("k")
    Ep_return = Ep["power", flip_idx:].C.set_plot_color("r")
    fl.plot(
        Ep_progressive,
        "o",
        label="progressive sat.",
        human_units=False,
    )
    fl.plot(
        Ep_return,
        "s",
        label="return check",
        human_units=False,
    )
    fl.plot(Ep_fit_curve, "k:", alpha=0.6, label="fit", human_units=False)
    ax = plt.gca()
    ax.text(
        0.5,
        0.70,
        "\n".join(
            [
                f"${sp.latex(sp.Symbol(j))} = {k:0.5g}$"
                for j, k in Ep_fit.output().items()
            ]
            + [r"$k_{\sigma} = %0.5g\ \mathrm{M^{-1}s^{-1}}$" % ksig]
        ),
        ha="center",
        va="center",
        size=9,
        transform=ax.transAxes,
    )
    plt.xlabel("Power / W")
    plt.ylabel(r"$E(p)$")
    plt.legend()
    # }}}

    # {{{ T₁(p) map
    fl.next(f"{sample_label} – T₁(p)")
    T1p = 1.0 / R1p
    fl.plot(T1p, "ko", label=r"$T_1(p)$")
    plt.xlabel("Power / W")
    plt.ylabel(r"$T_1$ / s")
    plt.legend()
    # }}}

    # {{{ R₁(p) with k_ρ⁻¹ linear fit
    fl.next(r"$R_1(p)$")
    fl.plot(R1p, "ko", label="experimental")
    fl.plot(
        R1p_fit,
        "k-",
        alpha=0.5,
        label=(
            r"linear $k_\rho^{-1}$ fit"
            + (" [exact—2 pts]" if n_fir == 2 else "")
        ),
    )
    plt.xlabel("Power / W")
    plt.ylabel(r"$R_1$ / s$^{-1}$")
    plt.legend()
    # }}}

    # {{{ k_ρ(p)
    fl.next(r"$k_\rho(p)$")
    fl.plot(krho, "ko", label=r"$k_\rho$ experimental")
    fl.plot(krho_fit, "k-", alpha=0.5, label=r"from linear $k_\rho^{-1}$ fit")
    plt.xlabel("Power / W")
    plt.ylabel(r"$k_\rho$ / M$^{-1}$ s$^{-1}$")
    plt.legend()
    # }}}

    # {{{ motional narrowing plots (only when enough FIR points)
    if do_mn_fit:
        fl.next(r"$R_{1,w}(p)$ – motional narrowing")
        fl.plot(R10p_for_mn, "ko", label=r"$R_{1,w}(p)$")
        fl.plot(R10p_mn_fit, "k-", alpha=0.5, label=r"Arrhenius $\tau_c$ fit")
        plt.xlabel("Power / W")
        plt.ylabel(r"$R_{1,w}$ / s$^{-1}$")
        plt.legend()

        fl.next(r"$k_\rho(p)$ – motional narrowing")
        fl.plot(krho, "ko", label=r"$k_\rho$ experimental")
        fl.plot(krho_mn_fit, "k-", alpha=0.5, label=r"Arrhenius $\tau_c$ fit")
        plt.xlabel("Power / W")
        plt.ylabel(r"$k_\rho$ / M$^{-1}$ s$^{-1}$")
        plt.legend()
    else:
        print(
            f"\nMotional narrowing plots skipped "
            f"(need ≥ {MIN_FIR_FOR_MN} FIR points, have {n_fir}). "
            f"Add more FIR_ nodes to the HDF5 file — no code changes needed."
        )
    # }}}
