r"""Fit ODNP data from saved integral tables.

Run ``generate_integrals_RealData.py`` first.  This script loads the
top-level ``Ep``, ``R1p``, and ``T1p`` nodes from the generated integrals
file and follows the ODNP book-chapter analysis:

.. math::

    \epsilon(p) &= 1-E(p) \\
    k_\rho(p) &= \frac{R_1(p)-R_{1,0}(p)}{C_\mathrm{SL}} \\
    k_\sigma s(p) &=
        \frac{\epsilon(p)R_1(p)}{C_\mathrm{SL}}
        \left|\frac{\omega_H}{\omega_e}\right|.

The measured :math:`k_\rho^{-1}(p)` values are fit to an
uncertainty-weighted low-order polynomial and reinserted into the
relaxation expression before fitting the saturation curve.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd
import sympy as sp

# {{{ changeable parameters
# thisfile = "260724_TMTPDI_ODNP_1.h5"
thisfile = "260724_TMTPDI_ODNP_1.h5"
output_dir = Path("/Users/atahan/exp_data/Atahan_Processed_Data/ODNP")
dataset_id = thisfile.removesuffix(".h5")
output_file = f"{dataset_id}_integrals.h5"
KRHO_INV_POLY_ORDER = 2
PHALF_MIN_W = 0.05
PHALF_MAX_POWER_FACTOR = 2.0
PLOT_FONT_SIZE = 18
# }}}


def main():
    # {{{ Load and validate the generated integral tables
    Ep = psd.nddata_hdf5(
        f"{output_file}/Ep",
        directory=str(output_dir),
    )
    R1p = psd.nddata_hdf5(
        f"{output_file}/R1p",
        directory=str(output_dir),
    )
    T1p = psd.nddata_hdf5(
        f"{output_file}/T1p",
        directory=str(output_dir),
    )
    for table_name, table in (("Ep", Ep), ("R1p", R1p), ("T1p", T1p)):
        if "power" not in table.dimlabels:
            raise ValueError(f"{table_name} does not have a power dimension")
        if not np.isfinite(table.data).all():
            raise ValueError(
                f"{table_name} contains non-finite values; regenerate "
                "the integral tables"
            )
    if not np.allclose(T1p.data, 1.0 / R1p.data, rtol=1e-7, atol=1e-10):
        raise ValueError("T1p is not the reciprocal of R1p")

    acq = Ep.get_prop("acq_params")
    if acq is None:
        raise ValueError("Ep is missing acq_params metadata")
    required_acq_params = (
        "concentration",
        "T1water_cold",
        "T1water_hot",
        "max_power",
        "guessed_MHz_to_GHz",
    )
    missing_acq_params = [k for k in required_acq_params if k not in acq]
    if missing_acq_params:
        raise KeyError(
            "acq_params is missing " + ", ".join(missing_acq_params)
        )
    concentration = float(acq["concentration"])
    if concentration <= 0:
        raise ValueError("concentration must be positive")
    chemical = acq.get("chemical", b"TEMPOL")
    if isinstance(chemical, bytes):
        chemical = chemical.decode()
    sample_label = f"{1e3 * concentration:g} mM {chemical}"
    # }}}

    # {{{ Determine progressive enhancement powers
    # Enhancement powers rise and then return to lower powers as a
    # reproducibility check.  The maximum-power point belongs to the
    # progressive series.
    flip_idx = np.argmax(Ep.getaxis("power")) + 1
    progressive_power = Ep["power", :flip_idx].getaxis("power")
    first_nonzero_power = np.flatnonzero(progressive_power > 0)
    if len(first_nonzero_power) == 0:
        raise ValueError(
            "No nonzero progressive powers are available for fitting"
        )
    first_nonzero_power = first_nonzero_power[0]
    # }}}

    # {{{ Generate water R1,0(p) and select physical R1(p) points
    T10_p = np.r_[
        acq["T1water_cold"],
        (acq["T1water_hot"] - acq["T1water_cold"]) / acq["max_power"],
    ]
    R10_p_all = 1.0 / R1p.fromaxis("power").eval_poly(
        T10_p,
        "power",
    )
    R1p_values = R1p.real.data
    R1p_outlier_mask = (
        ~np.isfinite(R1p_values)
        | ~np.isfinite(R10_p_all.real.data)
        | (R1p_values <= R10_p_all.real.data)
    )
    for bound_name, comparison in (
        ("R1_fit_lower_bound", np.less),
        ("R1_fit_upper_bound", np.greater),
    ):
        R1p_bound = R1p.get_prop(bound_name)
        if R1p_bound is not None:
            R1p_bound = float(R1p_bound)
            R1p_outlier_mask |= np.isclose(
                R1p_values,
                R1p_bound,
                rtol=1e-5,
                atol=1e-8,
            ) | comparison(R1p_values, R1p_bound)
    R1p_fit_idx = np.flatnonzero(~R1p_outlier_mask)
    R1p_outlier_idx = np.flatnonzero(R1p_outlier_mask)
    if len(R1p_fit_idx) <= KRHO_INV_POLY_ORDER:
        raise ValueError(
            "Not enough physical R1p points for the requested "
            f"order-{KRHO_INV_POLY_ORDER} k_rho inverse fit"
        )
    R1p_for_fit = R1p["power", R1p_fit_idx]
    R10_p = R10_p_all["power", R1p_fit_idx]
    if len(R1p_outlier_idx) > 0:
        R1p_outliers = R1p["power", R1p_outlier_idx]
        print(
            "Excluding R1p point(s) at powers "
            f"{R1p_outliers.getaxis('power')} W: "
            f"{R1p_outliers.data} s^-1"
        )
    else:
        R1p_outliers = None
    # }}}

    # {{{ Fit k_rho inverse and reconstruct R1(p)
    krho = (R1p_for_fit - R10_p) / concentration
    krho_inv = 1.0 / krho
    krho_inv_error = krho_inv.get_error()
    if krho_inv_error is None:
        raise ValueError(
            "R1p has no uncertainty data for weighted k_rho fitting"
        )
    finite_krho = (
        np.isfinite(krho_inv.real.data)
        & np.isfinite(krho_inv_error)
        & (krho_inv_error > 0)
    )
    if finite_krho.sum() <= KRHO_INV_POLY_ORDER:
        raise ValueError("Not enough finite weighted k_rho points for fitting")
    krho_inv_coeff = np.polynomial.polynomial.polyfit(
        R1p_for_fit.getaxis("power")[finite_krho],
        krho_inv.real.data[finite_krho],
        deg=KRHO_INV_POLY_ORDER,
        w=1.0 / krho_inv_error[finite_krho],
    )
    A, KsigmaSmax, phalf, power = sp.symbols(
        "A KsigmaSmax phalf power",
        real=True,
    )
    krho_inv_poly_expr = sum(
        krho_inv_coeff[j] * power**j for j in range(KRHO_INV_POLY_ORDER + 1)
    )
    R1p_expr = (T10_p[0] + T10_p[1] * power) ** -1 + (
        concentration / krho_inv_poly_expr
    )
    p_max = max(
        Ep.getaxis("power").max(),
        R1p.getaxis("power").max(),
    )
    phalf_max = PHALF_MAX_POWER_FACTOR * p_max
    powers_fine = psd.nddata(
        np.r_[0:p_max:300j],
        "power",
    ).set_units("power", "W")
    R1p_func = sp.lambdify(power, R1p_expr, "numpy")
    krho_inv_func = sp.lambdify(
        power,
        krho_inv_poly_expr,
        "numpy",
    )
    R1p_fit = powers_fine.fromaxis("power").run(R1p_func)
    krho_zero = float(1.0 / krho_inv_func(0.0))
    krho_hot = float(1.0 / krho_inv_func(p_max))
    # }}}

    # {{{ Fit k_sigma s(p)
    saturation_expr = power / (power + phalf)
    omegaH_over_omegaE = float(acq["guessed_MHz_to_GHz"]) * 1e-3
    epsilon = 1.0 - Ep
    R1_at_Ep = Ep.fromaxis("power").run(R1p_func)
    ksigma_s = epsilon * R1_at_Ep * omegaH_over_omegaE / concentration
    ksigma_s_for_fit = ksigma_s["power", first_nonzero_power:flip_idx]
    ksigma_s_fit = psd.lmfitdata(ksigma_s_for_fit.real)
    ksigma_s_fit.functional_form = KsigmaSmax * saturation_expr
    KsigmaSmax_guess = float(np.nanmax(ksigma_s_for_fit.real.data))
    if not np.isfinite(KsigmaSmax_guess) or KsigmaSmax_guess <= 0:
        KsigmaSmax_guess = abs(float(ksigma_s_for_fit.real.data[-1]))
    if KsigmaSmax_guess <= 0:
        raise ValueError("Cannot generate a positive k_sigma guess")
    ksigma_s_fit.set_guess(
        KsigmaSmax=dict(
            value=KsigmaSmax_guess,
            min=0.0,
            max=3.0 * KsigmaSmax_guess,
        ),
        phalf=dict(
            value=float(acq.get("guessed_phalf", 0.2)),
            min=PHALF_MIN_W,
            max=phalf_max,
        ),
    )
    # lmfitdata uses its symbolic Jacobian and attached data errors.
    ksigma_s_fit.fit()
    ksigma_s_fit_curve = ksigma_s_fit.eval(100)
    ksigma = float(ksigma_s_fit.output("KsigmaSmax"))
    phalf_value = float(ksigma_s_fit.output("phalf"))
    coupling_factor = ksigma / krho_zero
    # }}}

    # {{{ Fit normalized E(p) using the branch's legacy expression
    Ep_fit = psd.lmfitdata(Ep["power", :flip_idx].real)
    M0 = Ep["power", 0].real.item()
    Ep_fit.functional_form = M0 - (M0 * A * saturation_expr / R1p_expr)
    phalf_guess = float(acq.get("guessed_phalf", 0.2))
    progressive_Ep = Ep["power", :flip_idx].real
    max_epsilon_idx = np.nanargmax(1.0 - progressive_Ep.data)
    p_for_A_guess = progressive_power[max_epsilon_idx]
    s_for_A_guess = p_for_A_guess / (p_for_A_guess + phalf_guess)
    if s_for_A_guess <= 0:
        raise ValueError("Cannot generate A guess from zero power")
    A_guess = float(
        np.real(
            (1.0 - progressive_Ep.data[max_epsilon_idx] / M0)
            * R1p_func(p_for_A_guess)
            / s_for_A_guess
        )
    )
    Ep_fit.set_guess(
        A=dict(
            value=A_guess,
            min=0.0,
            max=3.0 * A_guess,
        ),
        phalf=dict(
            value=phalf_guess,
            min=PHALF_MIN_W,
            max=phalf_max,
        ),
    )
    Ep_fit.fit()
    Ep_fit_curve = Ep_fit.eval(100)
    # }}}

    # {{{ Plot ODNP fits
    plt.rcParams.update(
        {
            "font.size": PLOT_FONT_SIZE,
            "axes.titlesize": PLOT_FONT_SIZE,
            "axes.labelsize": PLOT_FONT_SIZE,
            "xtick.labelsize": PLOT_FONT_SIZE,
            "ytick.labelsize": PLOT_FONT_SIZE,
            "legend.fontsize": plt.rcParamsDefault["font.size"],
            "figure.titlesize": PLOT_FONT_SIZE,
        }
    )
    with psd.figlist_var() as fl:
        fl.basename = output_file
        fig = plt.figure(figsize=(10, 7.5), layout="constrained")
        fig.suptitle(sample_label)
        gs = fig.add_gridspec(2, 2)
        ax_epsilon = fig.add_subplot(gs[0, 0])
        ax_R1 = fig.add_subplot(gs[0, 1])
        ax_ksigma = fig.add_subplot(gs[1, :])
        fl.next("ODNP summary", fig=fig)

        psd.plot(
            epsilon["power", :flip_idx].C.set_plot_color("k"),
            "o",
            ax=ax_epsilon,
            label="progressive",
            human_units=False,
        )
        if flip_idx < Ep.shape["power"]:
            psd.plot(
                epsilon["power", flip_idx:].C.set_plot_color("r"),
                "s",
                ax=ax_epsilon,
                label="return check",
                human_units=False,
            )
        psd.plot(
            (1.0 - Ep_fit_curve).C.set_plot_color("k"),
            ":",
            ax=ax_epsilon,
            alpha=0.45,
            label="legacy raw E fit mapped to epsilon",
            human_units=False,
        )
        ax_epsilon.set_xlabel("Power / W")
        ax_epsilon.set_ylabel(r"$\epsilon(p) = 1 - E(p)$")
        ax_epsilon.legend()
        psd.gridandtick(ax_epsilon)

        psd.plot(
            R1p_for_fit.C.set_plot_color("k"),
            "o",
            ax=ax_R1,
            label="fit points",
            human_units=False,
        )
        psd.plot(
            R1p_fit.C.set_plot_color("k"),
            "-",
            ax=ax_R1,
            alpha=0.5,
            label=(
                rf"weighted order-{KRHO_INV_POLY_ORDER} " r"$k_\rho^{-1}$ fit"
            ),
            human_units=False,
        )
        R1_plot_data = np.r_[R1p_for_fit.real.data, R1p_fit.real.data]
        R1_y_pad = 0.1 * np.ptp(R1_plot_data)
        if R1_y_pad == 0:
            R1_y_pad = 1.0
        R1_y_limits = (
            np.nanmin(R1_plot_data) - R1_y_pad,
            np.nanmax(R1_plot_data) + R1_y_pad,
        )
        if R1p_outliers is not None:
            R1p_outlier_plot = R1p_outliers.C.set_plot_color("r")
            R1p_outlier_plot.data = np.clip(
                R1p_outlier_plot.data,
                R1_y_limits[0],
                R1_y_limits[1],
            )
            R1p_outlier_plot.set_error(None)
            psd.plot(
                R1p_outlier_plot,
                "s",
                ax=ax_R1,
                label="outlier, excluded from fit",
                human_units=False,
            )
        ax_R1.set_ylim(*R1_y_limits)
        ax_R1.set_xlabel("Power / W")
        ax_R1.set_ylabel(r"$R_1$ / s$^{-1}$")
        ax_R1.legend()
        psd.gridandtick(ax_R1)

        if first_nonzero_power > 0:
            psd.plot(
                ksigma_s["power", :first_nonzero_power].C.set_plot_color(
                    "0.5"
                ),
                "x",
                ax=ax_ksigma,
                human_units=False,
            )
        psd.plot(
            ksigma_s_for_fit.C.set_plot_color("k"),
            "o",
            ax=ax_ksigma,
            label=r"$k_\sigma s(p)$ fit points",
            human_units=False,
        )
        if flip_idx < Ep.shape["power"]:
            psd.plot(
                ksigma_s["power", flip_idx:].C.set_plot_color("r"),
                "s",
                ax=ax_ksigma,
                label="return check",
                human_units=False,
            )
        psd.plot(
            ksigma_s_fit_curve.C.set_plot_color("k"),
            "-",
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
                    (
                        r"$k_{\sigma}s_{\max} = %0.5g"
                        r"\ \mathrm{M^{-1}s^{-1}}$"
                    )
                    % ksigma,
                    r"$p_{1/2} = %0.5g\ \mathrm{W}$" % phalf_value,
                ]
            ),
            ha="right",
            va="bottom",
            size=PLOT_FONT_SIZE,
            transform=ax_ksigma.transAxes,
        )
        ax_ksigma.set_xlabel("Power / W")
        ax_ksigma.set_ylabel(r"$k_\sigma s(p)$ / M$^{-1}$ s$^{-1}$")
        ax_ksigma.legend()
        psd.gridandtick(ax_ksigma)
    # }}}

    # {{{ Console summary
    print(f"dataset: {thisfile}")
    print(f"sample: {sample_label}")
    print(f"pmax: {p_max:#0.6g} W")
    print(f"k_rho(0): {krho_zero:#0.6g} M^-1 s^-1")
    print(f"k_rho(pmax): {krho_hot:#0.6g} M^-1 s^-1")
    print(f"k_sigma: {ksigma:#0.6g} M^-1 s^-1")
    print(f"p_1/2: {phalf_value:#0.6g} W")
    print(f"coupling factor xi: {coupling_factor:#0.6g}")
    print(
        f"k_rho inverse polynomial coefficients: {np.asarray(krho_inv_coeff)}"
    )
    # }}}
    return {
        "krho_zero": krho_zero,
        "krho_hot": krho_hot,
        "ksigma": ksigma,
        "phalf": phalf_value,
        "pmax": p_max,
        "coupling_factor": coupling_factor,
        "krho_inv_coeff": np.asarray(krho_inv_coeff),
    }


if __name__ == "__main__":
    main()
