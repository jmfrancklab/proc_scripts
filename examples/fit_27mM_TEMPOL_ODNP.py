r"""High-concentration TEMPOL ODNP: T1 map and motional narrowing fits."""

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

thisfile, thisexptype, nodename = (
    "260526_hydroxytempo_ODNP_2.h5",
    "B27/ODNP",
    "ODNP",
)
C = 27e-3  # mol/L TEMPOL fallback
T0 = 298.15
heating_K_per_W = 5.0
Ea = 15e3
tau0 = 30e-12


def node_power(nodename):
    m = re.search(r"([0-9]+(?:p[0-9]+)?)dBm", nodename)
    if m:
        return prscr.dBm2power(float(m.group(1).replace("p", "."))).item()
    if "noPower" in nodename:
        return 0.0
    raise ValueError(f"Can't infer microwave power from node {nodename!r}")


def sorted_fir_nodes():
    filename = psd.search_filename(
        thisfile, exp_type=thisexptype, unique=True
    )
    with h5py.File(filename, "r") as f:
        nodes = [j for j in f.keys() if j.startswith("FIR_")]
    return sorted([(j, node_power(j)) for j in nodes], key=lambda x: x[1])


FIR_nodes = sorted_fir_nodes()

with psd.figlist_var() as fl:
    fl.basename = thisfile
    Ep = psd.find_file(
        thisfile,
        exp_type=thisexptype,
        expno=nodename,
        lookup=prscr.lookup_table,
        fl=fl,
    )
    if Ep.get_prop("log") is None:
        Ep = prscr.attach_log_data_from_file(Ep, thisfile, thisexptype)
    m = re.search(".*ODNP.*v([0-9]+)$", Ep.get_prop("postproc_type"))
    if m and int(m.groups()[0]) < 6:
        Ep = generate_coordinates_from_log(Ep, fl=fl)
    orig_axis = Ep["indirect"]
    orig_axis_error = Ep.get_error("indirect")
    Ep["indirect"] = orig_axis["time"]
    Ep.set_units("indirect", "s")
    Ep.set_error("indirect", orig_axis_error["time"])
    Ep, _ = prscr.rough_table_of_integrals(Ep, fl=fl)

Ep.set_error(Ep["indirect", 0].item() * 0.01)
Ep /= Ep["indirect", 0:1]
Ep.set_error("indirect", orig_axis_error["power"])
Ep["indirect"] = orig_axis["power"]
Ep.set_units("indirect", "W").rename("indirect", "power")

R1p = psd.ndshape([("power", len(FIR_nodes))]).alloc(dtype=np.float64)
R1p.setaxis("power", [j[1] for j in FIR_nodes]).set_units("power", "W")
acq = None

with psd.figlist_var() as fl:
    for j, (thisnodename, _) in enumerate(FIR_nodes):
        s = psd.find_file(
            thisfile,
            exp_type=thisexptype,
            expno=thisnodename,
            lookup=prscr.lookup_table,
        )
        if "nScans" in s.dimlabels:
            s = prscr.clock_correct(s)
        s = s.squeeze()
        acq = s.get_prop("acq_params")
        # C = acq.get("concentration", C)
        C = C
        s, ax_last = prscr.rough_table_of_integrals(
            s, fl=fl, signal_range=(-1500, 500)
        )
        Mi, R1, vd = sp.symbols("M_inf R_1 vd", real=True)
        W = (
            s.get_prop("acq_params")["FIR_rep_us"] * 1e-6
            + s.get_prop("acq_params")["acq_time_ms"] * 1e-3
        )
        s = psd.lmfitdata(s)
        s.functional_form = Mi * (1 - (2 - sp.exp(-W * R1)) * sp.exp(-vd * R1))
        scale = 10 ** psd.det_unit_prefactor(s.get_units("vd"))
        s.set_guess(
            M_inf=dict(
                value=s.max().item(),
                min=0.1 * s.max().item(),
                max=1.5 * s.max().item(),
            ),
            R_1=dict(value=0.8 * scale, min=0.01 * scale, max=100 * scale),
        )
        s.fit()
        R1p["power", j] = s.output("R_1")
        psd.plot(s.eval(200), ax=ax_last, alpha=0.5)

flip_idx = np.where(np.diff(Ep.getaxis("power")) < 0)[0][0] + 1
sample_label = f"{1e3 * C:g} mM {acq.get('chemical', 'TEMPOL')}"

T10_p = np.r_[
    acq["T1water_cold"],
    (acq["T1water_hot"] - acq["T1water_cold"]) / acq["max_power"],
]
R10_p = 1 / R1p.fromaxis("power").eval_poly(T10_p, "power")
T1p = 1 / R1p
krho = (R1p - R10_p) / C
krho_inv = C / (R1p - R10_p)
krho_inv_coeff = krho_inv.polyfit("power", order=1)
powers_fine = psd.nddata(
    np.r_[0 : max(Ep.getaxis("power").max(), R1p.getaxis("power").max()) : 300j],
    "power",
)

M0, A, phalf, p = sp.symbols("M0 A phalf power", real=True)
R1p_expr = (T10_p[0] + T10_p[1] * p) ** -1 + (
    C / (krho_inv_coeff[0] + krho_inv_coeff[1] * p)
)
R1p_fit = psd.nddata(sp.lambdify(p, R1p_expr)(powers_fine.data), "power").setaxis(
    "power", powers_fine.data
)
krho_fit = psd.nddata(
    1 / (krho_inv_coeff[0] + krho_inv_coeff[1] * powers_fine.data),
    "power",
).setaxis("power", powers_fine.data)

mn_tau_expr = tau0 * sp.exp(Ea / 8.314 * (1 / (T0 + heating_K_per_W * p) - 1 / T0))
offset, scale = sp.symbols("offset scale", real=True)
mn_form = offset + scale * mn_tau_expr


def fit_motional_narrowing(y, p_axis, guess_scale):
    f = psd.lmfitdata(y.real)
    f.functional_form = mn_form
    f.set_guess(offset=y.real.min().item(), scale=guess_scale)
    f.fit()
    return f.eval(p_axis)


T10p_for_mn = R1p.fromaxis("power").eval_poly(T10_p, "power")
R10p_for_mn = 1 / T10p_for_mn
R10p_mn_fit = fit_motional_narrowing(R10p_for_mn, powers_fine.data, 1e10)
krho_mn_fit = fit_motional_narrowing(krho, powers_fine.data, 1e13)

sp_expr = p / (p + phalf)
Ep_fit = psd.lmfitdata(Ep["power", :flip_idx].real)
Ep_fit.functional_form = M0 - ((M0 * A * sp_expr) / R1p_expr)
A_guess = (
    1
    - (
        R1p["power", 0].real.item()
        * (Ep["power", flip_idx].real.item() / Ep["power", 0].real.item())
    ).real
)
Ep_fit.set_guess(
    M0=dict(value=Ep["power", 0].real.item(), min=0.1, max=2.0),
    A=dict(value=A_guess, min=0.2 * A_guess, max=3 * A_guess),
    phalf=dict(value=acq.get("guessed_phalf", 0.2), min=0.05, max=1.0),
)
Ep_fit.fit()
Ep_fit_curve = Ep_fit.eval(100)
ksig = Ep_fit.output("A") * acq["guessed_MHz_to_GHz"] * 1e-3 / C

with psd.figlist_var() as fl:
    fl.next(f"{sample_label} enhancement")
    fl.skip_units_check()
    fl.plot(Ep["power", :flip_idx], "ko", label="up")
    fl.plot(Ep["power", flip_idx:], "ro", label="return")
    fl.plot(Ep_fit_curve, "k:", alpha=0.5, label="fit")
    ax = plt.gca()
    ax.text(
        0.5,
        0.7,
        "\n".join(
            [f"${sp.latex(sp.Symbol(j))} = {k:0.5g}$" for j, k in Ep_fit.output().items()]
            + [r"$k_{\sigma} = %0.6f M^{-1}s^{-1}$" % ksig]
        ),
        ha="center",
        va="center",
        size=10,
        transform=ax.transAxes,
    )
    plt.xlabel("Power / W")
    plt.ylabel(r"$M_0 E(p)$")

    fl.next(f"{sample_label} T1 map")
    fl.plot(T1p, "o", label=r"$T_1(p)$")
    plt.xlabel("Power / W")
    plt.ylabel(r"$T_1$ / s")

    fl.next(r"$R_1(p)$ from $k_\rho^{-1}$ line")
    fl.skip_units_check()
    fl.plot(R1p, "o", label="experimental")
    fl.plot(R1p_fit, "k-", alpha=0.5, label="fit")
    plt.ylabel(r"$R_1$ / s$^{-1}$")

    fl.next(r"$k_\rho$")
    fl.skip_units_check()
    fl.plot(krho, "o", label=r"$k_\rho$")
    fl.plot(krho_fit, "k-", alpha=0.5, label=r"from linear $k_\rho^{-1}$")
    plt.ylabel(r"$k_\rho$ / M$^{-1}$ s$^{-1}$")

    fl.next("water motional narrowing")
    fl.skip_units_check()
    fl.plot(R10p_for_mn, "o", label=r"$R_{1,w}$")
    fl.plot(R10p_mn_fit, "k-", alpha=0.5, label="fit")
    plt.ylabel(r"$R_{1,w}$ / s$^{-1}$")

    fl.next(r"$k_\rho$ motional narrowing")
    fl.skip_units_check()
    fl.plot(krho, "o", label=r"$k_\rho$")
    fl.plot(krho_mn_fit, "k-", alpha=0.5, label="fit")
    plt.ylabel(r"$k_\rho$ / M$^{-1}$ s$^{-1}$")
