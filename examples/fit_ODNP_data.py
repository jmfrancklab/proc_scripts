r"""Fit ODNP data from saved integral tables.

Run ``generate_integrals_RealData.py`` first.  This script loads the saved
``Ep`` and ``R1p`` tables and applies the same ODNP model used in
``fit_27mM_TEMPOL_ODNP.py``:

    E(p) = M0 - M0 A s(p) / R1(p)

with ``R1(p)`` built from the water T1 interpolation and a linear fit to
``k_rho^{-1}(p)``.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd
import sympy as sp


# {{{ user-editable block
output_dir = Path(__file__).resolve().parent / "generated"
output_file = "ODNP_integrals_27mM.h5"
dataset_id = "260615_hydroxytempo_ODNP_3"

T0 = 298.15
heating_K_per_W = 5.0
Ea = 15e3
tau0 = 30e-12
MIN_FIR_FOR_MN = 4
# }}}


def load_table(node_name):
    h5path = f"{output_file}/{dataset_id}/{node_name}"
    if not (output_dir / output_file).exists():
        raise FileNotFoundError(
            f"Missing {output_dir / output_file}. "
            "Run examples/generate_integrals_RealData.py first."
        )
    return psd.nddata_hdf5(h5path, directory=str(output_dir))


def fit_motional_narrowing(y_nddata, p_fine, mn_form, guess_scale):
    f = psd.lmfitdata(y_nddata.real)
    f.functional_form = mn_form
    f.set_guess(
        offset=dict(value=float(y_nddata.real.data.min()), vary=True),
        scale=dict(value=guess_scale, vary=True),
    )
    f.fit()
    if hasattr(p_fine, "data"):
        p_fine = p_fine.data
    return f.eval(p_fine)


Ep = load_table("Ep")
R1p = load_table("R1p")
if np.isnan(R1p.data).any():
    raise ValueError("R1p table contains NaNs; regenerate integrals first.")

acq = Ep.get_prop("acq_params")
C = float(acq["concentration"])
chemical = acq.get("chemical", b"TEMPOL")
if isinstance(chemical, bytes):
    chemical = chemical.decode()
sample_label = f"{1e3 * C:g} mM {chemical}"

flip_idx = np.argmax(Ep.getaxis("power")) + 1
T10_p = np.r_[
    acq["T1water_cold"],
    (acq["T1water_hot"] - acq["T1water_cold"]) / acq["max_power"],
]
R10_p = 1.0 / R1p.fromaxis("power").eval_poly(T10_p, "power")
krho = (R1p - R10_p) / C
krho_inv = C / (R1p - R10_p)
krho_inv_coeff = krho_inv.polyfit("power", order=1)

M0, A, phalf, p = sp.symbols("M0 A phalf power", real=True)
R1p_expr = (T10_p[0] + T10_p[1] * p) ** -1 + (
    C / (krho_inv_coeff[0] + krho_inv_coeff[1] * p)
)
p_max = max(Ep.getaxis("power").max(), R1p.getaxis("power").max())
powers_fine = psd.nddata(np.r_[0:p_max:300j], "power").set_units("power", "W")
R1p_fit = powers_fine.fromaxis("power").run(sp.lambdify(p, R1p_expr, "numpy"))
krho_fit = powers_fine.fromaxis("power").run(
    lambda x: 1.0 / (krho_inv_coeff[0] + krho_inv_coeff[1] * x)
)

do_mn_fit = len(R1p.getaxis("power")) >= MIN_FIR_FOR_MN
if do_mn_fit:
    offset_sym, scale_sym = sp.symbols("offset scale", real=True)
    mn_tau_expr = tau0 * sp.exp(
        Ea / 8.314 * (1.0 / (T0 + heating_K_per_W * p) - 1.0 / T0)
    )
    mn_form = offset_sym + scale_sym * mn_tau_expr
    R10p_for_mn = 1.0 / R1p.fromaxis("power").eval_poly(T10_p, "power")
    R10p_mn_fit = fit_motional_narrowing(
        R10p_for_mn,
        powers_fine,
        mn_form,
        1e10,
    )
    krho_mn_fit = fit_motional_narrowing(krho, powers_fine, mn_form, 1e13)

sp_expr = p / (p + phalf)
Ep_fit = psd.lmfitdata(Ep["power", :flip_idx].real)
Ep_fit.functional_form = M0 - (M0 * A * sp_expr) / R1p_expr
A_guess = 1.0 - (
    R1p["power", 0].real.item()
    * (Ep["power", flip_idx].real.item() / Ep["power", 0].real.item())
)
M0_val = Ep["power", 0].real.item()
Ep_fit.set_guess(
    M0=dict(value=M0_val, min=0.1 * M0_val, max=2.0 * M0_val),
    A=dict(value=A_guess, min=0.2 * A_guess, max=3.0 * A_guess),
    phalf=dict(value=float(acq.get("guessed_phalf", 0.2)), min=0.05, max=1.0),
)
Ep_fit.fit()
Ep_fit_curve = Ep_fit.eval(100)
ksig = Ep_fit.output("A") * float(acq["guessed_MHz_to_GHz"]) * 1e-3 / C

with psd.figlist_var() as fl:
    fl.basename = output_file

    fl.next(f"{sample_label} - E(p)")
    fl.plot(
        Ep["power", :flip_idx].C.set_plot_color("k"),
        "o",
        label="progressive",
        human_units=False,
    )
    fl.plot(
        Ep["power", flip_idx:].C.set_plot_color("r"),
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

    fl.next(f"{sample_label} - T1(p)")
    fl.plot(1.0 / R1p, "ko", label=r"$T_1(p)$", human_units=False)
    plt.xlabel("Power / W")
    plt.ylabel(r"$T_1$ / s")
    plt.legend()

    fl.next(r"$R_1(p)$")
    fl.plot(R1p, "ko", label="experimental", human_units=False)
    fl.plot(R1p_fit, "k-", alpha=0.5, label=r"linear $k_\rho^{-1}$ fit")
    plt.xlabel("Power / W")
    plt.ylabel(r"$R_1$ / s$^{-1}$")
    plt.legend()

    fl.next(r"$k_\rho(p)$")
    fl.plot(krho, "ko", label=r"$k_\rho$ experimental", human_units=False)
    fl.plot(
        krho_fit,
        "k-",
        alpha=0.5,
        label=r"from linear $k_\rho^{-1}$ fit",
        human_units=False,
    )
    plt.xlabel("Power / W")
    plt.ylabel(r"$k_\rho$ / M$^{-1}$ s$^{-1}$")
    plt.legend()

    if do_mn_fit:
        fl.next(r"$R_{1,w}(p)$ - motional narrowing")
        fl.plot(R10p_for_mn, "ko", label=r"$R_{1,w}(p)$", human_units=False)
        fl.plot(
            R10p_mn_fit,
            "k-",
            alpha=0.5,
            label=r"Arrhenius $\tau_c$ fit",
            human_units=False,
        )
        plt.xlabel("Power / W")
        plt.ylabel(r"$R_{1,w}$ / s$^{-1}$")
        plt.legend()

        fl.next(r"$k_\rho(p)$ - motional narrowing")
        fl.plot(krho, "ko", label=r"$k_\rho$ experimental", human_units=False)
        fl.plot(
            krho_mn_fit,
            "k-",
            alpha=0.5,
            label=r"Arrhenius $\tau_c$ fit",
            human_units=False,
        )
        plt.xlabel("Power / W")
        plt.ylabel(r"$k_\rho$ / M$^{-1}$ s$^{-1}$")
        plt.legend()
