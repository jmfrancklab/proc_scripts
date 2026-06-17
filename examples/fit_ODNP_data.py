r"""Fit ODNP data from saved integral tables.

Run ``generate_integrals_RealData.py`` first.  This script loads the top-level
``Ep``, ``R1p``, and ``T1p`` nodes from ``<source filename>_integrals.h5`` using
``find_file``, then follows the ODNP analysis in the DCCT paper and ODNP book
chapter:

    epsilon(p) = 1 - E(p)
    k_sigma s(p) = epsilon(p) R1(p) / C_SL * |omega_H / omega_e|
    k_sigma s(p) = k_sigma s_max p / (p_1/2 + p)

``R1(p)`` is interpolated by fitting ``k_rho^{-1}(p)`` to a low-order
polynomial, then reinserting that result into the relaxation expression,
matching the ODNP book chapter protocol.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd
import sympy as sp


# {{{ user-editable block
thisfile = "260615_hydroxytempo_ODNP_3.h5"
output_dir = Path("/Users/atahan/exp_data/Atahan_Processed_Data/ODNP")
dataset_id = thisfile.removesuffix(".h5")
output_file = f"{dataset_id}_integrals.h5"

KRHO_INV_POLY_ORDER = 2
# }}}


def load_table(node_name):
    """Load one top-level generated table using the broken example's style."""
    return psd.find_file(
        output_file,
        exp_type=str(output_dir),
        expno=node_name,
    )


# {{{ Load saved top-level integral tables using find_file
Ep = load_table("Ep")
R1p = load_table("R1p")
T1p = load_table("T1p")

if np.isnan(R1p.data).any():
    raise ValueError("R1p table contains NaNs; regenerate integrals first.")

R1p_values = R1p.real.data
R1p_upper_bound = R1p.get_prop("R1_fit_upper_bound")
if R1p_upper_bound is None:
    raise ValueError(
        "R1p table is missing R1_fit_upper_bound metadata; regenerate "
        "integrals with examples/generate_integrals_RealData.py."
    )
R1p_upper_bound = float(R1p_upper_bound)
R1p_outlier_mask = np.isclose(
    R1p_values, R1p_upper_bound, rtol=1e-8, atol=1e-8
) | (R1p_values > R1p_upper_bound)
R1p_fit_idx = np.flatnonzero(~R1p_outlier_mask)
R1p_outlier_idx = np.flatnonzero(R1p_outlier_mask)
if len(R1p_fit_idx) < 2:
    raise ValueError("Need at least two non-outlier R1p points for fitting.")
R1p_for_fit = R1p["power", R1p_fit_idx]
if len(R1p_outlier_idx) > 0:
    R1p_outliers = R1p["power", R1p_outlier_idx]
    print(
        "Excluding R1p outlier(s) at powers "
        f"{R1p_outliers.getaxis('power')} W: {R1p_outliers.data} s^-1"
    )
else:
    R1p_outliers = None

acq = Ep.get_prop("acq_params")
C = float(acq["concentration"])
chemical = acq.get("chemical", b"TEMPOL")
if isinstance(chemical, bytes):
    chemical = chemical.decode()
sample_label = f"{1e3 * C:g} mM {chemical}"
# }}}

# {{{ The powers go up, then return for a reproducibility check
flip_idx = np.argmax(Ep.getaxis("power")) + 1
# }}}

# {{{ Generate water R1(p) from the cold/hot water T1 interpolation
T10_p = np.r_[
    acq["T1water_cold"],
    (acq["T1water_hot"] - acq["T1water_cold"]) / acq["max_power"],
]
R10_p = 1.0 / R1p_for_fit.fromaxis("power").eval_poly(T10_p, "power")
if R1p_outliers is not None:
    R10_p_outliers = 1.0 / R1p_outliers.fromaxis("power").eval_poly(
        T10_p, "power"
    )
# }}}

# {{{ Fit kᵨ⁻¹ to the book-chapter polynomial and apply to R1(p)
# Section 7.4 of the ODNP book chapter prescribes measuring R1(p) at a limited
# set of microwave powers, converting those points to k_rho^{-1}(p), fitting a
# first- or second-order polynomial, and reinserting the interpolated k_rho(p)
# into the R1(p) expression used for enhancement analysis.
krho = (R1p_for_fit - R10_p) / C
krho_inv = 1 / krho
if R1p_outliers is not None:
    krho_outliers = (R1p_outliers - R10_p_outliers) / C
    krho_inv_outliers = C / (R1p_outliers - R10_p_outliers)
krho_inv_coeff = krho_inv.polyfit("power", order=KRHO_INV_POLY_ORDER)

M0, A, KsigmaSmax, phalf, p = sp.symbols(
    "M0 A KsigmaSmax phalf power",
    real=True,
)
krho_inv_poly_expr = sum(
    krho_inv_coeff[j] * p**j for j in range(KRHO_INV_POLY_ORDER + 1)
)
R1p_expr = (T10_p[0] + T10_p[1] * p) ** -1 + (C / krho_inv_poly_expr)
p_max = max(Ep.getaxis("power").max(), R1p.getaxis("power").max())
powers_fine = psd.nddata(np.r_[0:p_max:300j], "power").set_units("power", "W")
R1p_func = sp.lambdify(p, R1p_expr, "numpy")
krho_inv_func = sp.lambdify(p, krho_inv_poly_expr, "numpy")
R1p_fit = powers_fine.fromaxis("power").run(R1p_func)
krho_inv_fit = powers_fine.fromaxis("power").run(krho_inv_func)
krho_fit = powers_fine.fromaxis("power").run(lambda x: 1.0 / krho_inv_func(x))
# }}}

# {{{ Fit k_sigma s(p)
saturation_expr = p / (p + phalf)
omegaH_over_omegaE = float(acq["guessed_MHz_to_GHz"]) * 1e-3
epsilon = 1.0 - Ep
R1_at_Ep = Ep.fromaxis("power").run(R1p_func)
ksigma_s = epsilon * R1_at_Ep * omegaH_over_omegaE / C

first_nonzero_power = np.flatnonzero(
    Ep["power", :flip_idx].getaxis("power") > 0
)
if len(first_nonzero_power) == 0:
    raise ValueError(
        "No nonzero progressive powers available for k_sigma s(p) fit"
    )
first_nonzero_power = first_nonzero_power[0]
ksigma_s_for_fit = ksigma_s["power", first_nonzero_power:flip_idx]
ksigma_s_fit = psd.lmfitdata(ksigma_s_for_fit.real)
ksigma_s_fit.functional_form = KsigmaSmax * saturation_expr
KsigmaSmax_guess = float(np.nanmax(ksigma_s_for_fit.real.data))
if not np.isfinite(KsigmaSmax_guess) or KsigmaSmax_guess <= 0:
    KsigmaSmax_guess = abs(float(ksigma_s_for_fit.real.data[-1]))
ksigma_s_fit.set_guess(
    KsigmaSmax=dict(
        value=KsigmaSmax_guess,
        min=0.0,
        max=3.0 * KsigmaSmax_guess,
    ),
    phalf=dict(value=float(acq.get("guessed_phalf", 0.2)), min=0.05, max=1.0),
)
ksigma_s_fit.fit()
ksigma_s_fit_curve = ksigma_s_fit.eval(100)
ksig = ksigma_s_fit.output("KsigmaSmax")
# }}}

# {{{ Fit raw E(p) integrals with the legacy expression for comparison
# Legacy/raw-integral fit used by examples/broken/fit_ODNP_data.py.  This is
# algebraically related, but the paper/book protocol above fits k_sigma s(p).
Ep_fit = psd.lmfitdata(Ep["power", :flip_idx].real)
Ep_fit.functional_form = M0 - (M0 * A * saturation_expr) / R1p_expr
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
# }}}

with psd.figlist_var() as fl:
    fl.basename = output_file

    # {{{ Plot epsilon(p), R1(p), and k_sigma s(p)
    fig = plt.figure(figsize=(10, 7.5), constrained_layout=True)
    fig.suptitle(sample_label)
    gs = fig.add_gridspec(2, 2)
    ax_eps = fig.add_subplot(gs[0, 0])
    ax_R1 = fig.add_subplot(gs[0, 1])
    ax_ksigma = fig.add_subplot(gs[1, :])
    fl.next("ODNP summary", fig=fig)

    psd.plot(
        epsilon["power", :flip_idx].C.set_plot_color("k"),
        "o",
        ax=ax_eps,
        label="progressive",
        human_units=False,
    )
    psd.plot(
        epsilon["power", flip_idx:].C.set_plot_color("r"),
        "s",
        ax=ax_eps,
        label="return check",
        human_units=False,
    )
    psd.plot(
        1.0 - Ep_fit_curve,
        "k:",
        ax=ax_eps,
        alpha=0.45,
        label="legacy raw E fit mapped to epsilon",
        human_units=False,
    )
    ax_eps.set_xlabel("Power / W")
    ax_eps.set_ylabel(r"$\epsilon(p) = 1 - E(p)$")
    ax_eps.legend()
    psd.gridandtick(ax_eps)

    psd.plot(
        R1p_for_fit,
        "ko",
        ax=ax_R1,
        label="fit points",
        human_units=False,
    )
    if R1p_outliers is not None:
        psd.plot(
            R1p_outliers.C.set_plot_color("r"),
            "s",
            ax=ax_R1,
            label="outlier, excluded from fit",
            human_units=False,
        )
    psd.plot(
        R1p_fit,
        "k-",
        ax=ax_R1,
        alpha=0.5,
        label=rf"order-{KRHO_INV_POLY_ORDER} $k_\rho^{{-1}}$ fit",
        human_units=False,
    )
    ax_R1.set_xlabel("Power / W")
    ax_R1.set_ylabel(r"$R_1$ / s$^{-1}$")
    ax_R1.legend()
    psd.gridandtick(ax_R1)

    if first_nonzero_power > 0:
        psd.plot(
            ksigma_s["power", :first_nonzero_power].C.set_plot_color("0.5"),
            "x",
            ax=ax_ksigma,
            human_units=False,
        )
    psd.plot(
        ksigma_s["power", first_nonzero_power:flip_idx].C.set_plot_color("k"),
        "ko",
        ax=ax_ksigma,
        label=r"$k_\sigma s(p)$ fit points",
        human_units=False,
    )
    psd.plot(
        ksigma_s["power", flip_idx:].C.set_plot_color("r"),
        "s",
        ax=ax_ksigma,
        label="return check",
        human_units=False,
    )
    psd.plot(
        ksigma_s_fit_curve,
        "k-",
        ax=ax_ksigma,
        alpha=0.6,
        label=r"$k_\sigma s_{\max}p/(p_{1/2}+p)$",
        human_units=False,
    )
    ax_ksigma.text(
        0.97,
        0.05,
        "\n".join(
            [
                r"$k_{\sigma}s_{\max} = %0.5g\ \mathrm{M^{-1}s^{-1}}$" % ksig,
                r"$p_{1/2} = %0.5g\ \mathrm{W}$"
                % ksigma_s_fit.output("phalf"),
            ]
        ),
        ha="right",
        va="bottom",
        size=9,
        transform=ax_ksigma.transAxes,
    )
    ax_ksigma.set_xlabel("Power / W")
    ax_ksigma.set_ylabel(r"$k_\sigma s(p)$ / M$^{-1}$ s$^{-1}$")
    ax_ksigma.legend()
    psd.gridandtick(ax_ksigma)
    # }}}
