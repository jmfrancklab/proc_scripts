"""
Process a field-swept ODNP enhancement experiment
=================================================

Process a field sweep acquired by ``run_field_sweep.py``.  The measured
Hall-probe readback is used as the experimental field coordinate, each NMR
echo is reduced to an integral with ``rough_table_of_integrals``, and the
resulting enhancement profile is normalized to the off-resonance sweep edges.

The acquisition is assumed to have been set up so that the proton carrier is
on resonance at the magnetic field corresponding to the measured microwave
resonance.  This common resonance field defines the experimental
NMR/EPR-frequency ratio and the field-to-EPR-frequency conversion.  Fits are
performed internally in gauss so that centers and linewidths can be compared
directly with CW-EPR fits.  The final plot is shown versus EPR offset in MHz:

    epr_offset = nu_EPR(B) - nu_MW.

The selected fit model can be changed below:

    ``"single_voigt"``
        One pseudo-Voigt line.  This is the general-purpose default.

    ``"three_voigt"``
        Three pseudo-Voigt lines.  Intended for resolved TEMPOL-like spectra.

    ``"polynomial"``
        Empirical smoothing only.  Useful when no physical line shape is
        appropriate, but its coefficients and linewidth are not physical.

    ``None``
        Process and plot without fitting.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp


plt.rcParams["image.aspect"] = "auto"
# sphinx_gallery_thumbnail_number = 2
plt.rcParams.update(
    {
        "errorbar.capsize": 2,
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.figsize": (6, 4),
    }
)


# {{{ changeable parameters
if "thisfile" not in globals():
    thisfile, exp_type, nodename = (
        "260713_hydroxytempo_field_sweep.h5",
        "B27/field_dependent",
        "field_sweep_1",
    )

# Choose one of "single_voigt", "three_voigt", "polynomial", or None.
fit_model = "single_voigt"

# Number of points from each sweep edge used for the off-resonance reference.
n_edge = 2

# Used only when fit_model == "polynomial".
poly_order = 8

# Initial line positions for the optional three-Voigt model.  These are from
# the independent 27 mM TEMPOL CW-EPR fit.  Only the line separations are used:
# the triplet is shifted so that the middle line is near the common resonance
# field calculated from the acquisition parameters.
three_voigt_centers_G = np.array([3381.40, 3398.23, 3415.16])
three_voigt_sigma_G = np.array([0.576, 0.582, 0.580])
three_voigt_lambda_G = np.array([3.396, 3.356, 3.481]) / 2
three_voigt_eta = np.array([0.93, 0.92, 0.93])
# }}}


if (
    "SPHINX_GALLERY_RUNNING" in os.environ
    and os.environ["SPHINX_GALLERY_RUNNING"] == "True"
):
    thisfile, exp_type, nodename = (
        "260713_field_sweep.h5",
        "B27/field_dependent",
        "field_sweep_1",
    )


def _require_acq_params(acq_params):
    "Check that the acquisition parameters needed for the calibration exist."
    if acq_params is None:
        raise ValueError("The dataset does not contain acq_params")
    for name in ("carrierFreq_MHz", "gamma_eff_MHz_G", "uw_dip_center_GHz"):
        if name not in acq_params:
            raise KeyError(f"acq_params does not contain {name}")


def _field_readback_G(s):
    """
    Return the measured Hall-probe readback in gauss.

    Newer field-sweep acquisitions store this in the ``field_readback_G``
    property.  Some older files keep field information inside the structured
    indirect axis, so this function accepts those names as a fallback.
    """
    field_readback = s.get_prop("field_readback_G")
    if field_readback is not None:
        field_G = np.asarray(field_readback["NUMPY_DATA"], dtype=float)
    else:
        indirect_axis = s["indirect"]
        field_name = next(
            (
                name
                for name in ("Field_G", "field_G", "Field", "field", "B")
                if getattr(indirect_axis, "dtype", None) is not None
                and indirect_axis.dtype.names is not None
                and name in indirect_axis.dtype.names
            ),
            None,
        )
        if field_name is None:
            raise ValueError(
                "Could not find field_readback_G or a field column in "
                "the indirect axis"
            )
        field_G = np.asarray(indirect_axis[field_name], dtype=float)

    # Some acquisition files historically stored kG values under a G-like
    # field name.  A median field below 100 is not physically plausible for
    # these X-band ODNP experiments, so interpret it as kG.
    if np.nanmedian(abs(field_G)) < 100:
        field_G = field_G * 1e3
    return field_G


def _calibration_from_acq(acq_params):
    """
    Calculate the common resonance field and frequency conversions.

    The experimentally determined frequency ratio is stored as MHz/GHz, which
    is numerically equal to parts per thousand.  The EPR field conversion is
    the effective conversion needed for this experiment, not a hard-coded
    free-electron constant.
    """
    resonance_field_G = (
        acq_params["carrierFreq_MHz"] / acq_params["gamma_eff_MHz_G"]
    )
    nmr_epr_freq_ratio = (
        acq_params["carrierFreq_MHz"] / acq_params["uw_dip_center_GHz"]
    )
    epr_field_conversion_MHz_G = (
        acq_params["uw_dip_center_GHz"] * 1e3 / resonance_field_G
    )
    return dict(
        resonance_field_G=resonance_field_G,
        nmr_epr_freq_ratio=nmr_epr_freq_ratio,
        epr_field_conversion_MHz_G=epr_field_conversion_MHz_G,
    )


def _pseudo_voigt(B, Bcenter, sigma, lambda_L, eta):
    "Absorption pseudo-Voigt line shape used by the lmfitdata models."
    return eta * lambda_L**2 / ((B - Bcenter) ** 2 + lambda_L**2) + (
        1 - eta
    ) * sp.exp(-((B - Bcenter) ** 2) / (2 * sigma**2))


def _pseudo_voigt_derivative(B, Bcenter, sigma, lambda_L, eta):
    "Derivative pseudo-Voigt line shape used for field-swept enhancement."
    return eta * (
        -2
        * lambda_L**2
        * (B - Bcenter)
        / ((B - Bcenter) ** 2 + lambda_L**2) ** 2
    ) + (1 - eta) * (
        -(B - Bcenter)
        / sigma**2
        * sp.exp(-((B - Bcenter) ** 2) / (2 * sigma**2))
    )


def _fwhm_pseudo_voigt(sigma, lambda_L):
    "Return the Olivero-Longbothum Voigt FWHM approximation in gauss."
    lorentzian_fwhm = 2 * lambda_L
    gaussian_fwhm = 2.35482 * sigma
    return 0.5346 * lorentzian_fwhm + np.sqrt(
        0.2166 * lorentzian_fwhm**2 + gaussian_fwhm**2
    )


def _edge_reference(s, n_edge):
    "Return edge indices and edge-signal diagnostics for normalization."
    if s.shape["B"] < 2 * n_edge + 3:
        raise ValueError(
            "The sweep does not contain enough points for edge normalization"
        )
    edge_idx = np.r_[0:n_edge, s.shape["B"] - n_edge : s.shape["B"]]
    left_edge_signal = s.data.real[:n_edge].mean()
    right_edge_signal = s.data.real[-n_edge:].mean()
    edge_signal = s.data.real[edge_idx].mean()
    if np.isclose(edge_signal, 0):
        raise ValueError(
            "Cannot normalize the field sweep because the edge signal is zero"
        )
    return edge_idx, edge_signal, left_edge_signal, right_edge_signal


def _largest_deviation(fit, baseline, dim="B"):
    "Find the fitted point with the largest deviation from baseline."
    position = fit.C.run(lambda x: abs(x - baseline)).argmax(dim).item()
    enhancement = fit[dim:position].item().real
    return position, enhancement


def _to_epr_offset_axis(data, conversion):
    "Copy a B-axis dataset and convert the copied axis to EPR offset in MHz."
    data = data.copy()
    data["B"] = (
        np.asarray(data["B"], dtype=float) - conversion["resonance_field_G"]
    ) * conversion["epr_field_conversion_MHz_G"]
    data.rename("B", "epr_offset")
    data.set_units("epr_offset", "MHz")
    return data


def fit_polynomial(s, order=6):
    """
    Empirically smooth the enhancement profile with a polynomial.

    This is not a physical EPR model and does not report a linewidth.  It is
    included as a robust fallback for irregular spectra where a Voigt model is
    visibly inappropriate.
    """
    if order < 2:
        raise ValueError("poly_order must be at least 2")
    if order >= s.shape["B"] - 2:
        raise ValueError(
            "poly_order is too high for the number of field-sweep points"
        )
    c_poly = s.polyfit("B", order)
    fit = s.eval_poly(c_poly, "B", npts=800)
    peak_field_G, peak_enhancement = _largest_deviation(fit, 1.0)
    return fit, dict(
        model="polynomial",
        baseline=1.0,
        center_field_G=None,
        linewidth_G=None,
        peak_field_G=peak_field_G,
        peak_enhancement=peak_enhancement,
    )


def fit_single_voigt(s):
    """
    Fit a single derivative pseudo-Voigt in the measured field coordinate.

    This follows the same fitting strategy as the earlier field-sweep script
    inspired by ``spin_trapping.py``: make a physically reasonable guess from
    the positive and negative lobes, fit only the center first, and only then
    release the rest of the parameters.  That staging is intentionally
    conservative and is usually more stable than optimizing the full line
    shape immediately.
    """
    B, E0, A, Bcenter, sigma, lambda_L, eta = sp.symbols(
        "B E_0 A Bcenter sigma lambda_L eta",
        real=True,
    )
    fitdata = psd.lmfitdata(s)
    fitdata.functional_form = E0 + A * _pseudo_voigt_derivative(
        B, Bcenter, sigma, lambda_L, eta
    )

    B_axis_G = np.asarray(s["B"], dtype=float)
    sweep_width_G = np.ptp(B_axis_G)
    enhancement_span = max(np.ptp(s.data.real), 0.1)
    baseline_guess = np.mean(
        np.r_[s.data.real[:n_edge], s.data.real[-n_edge:]]
    )
    y_from_baseline = s.data.real - baseline_guess
    max_idx = np.argmax(y_from_baseline)
    min_idx = np.argmin(y_from_baseline)
    linewidth_guess_G = max(
        abs(B_axis_G[min_idx] - B_axis_G[max_idx]) / 2,
        sweep_width_G / 8,
        0.1,
    )
    center_guess_G = 0.5 * (B_axis_G[max_idx] + B_axis_G[min_idx])
    amp_guess = (
        y_from_baseline[max_idx] - y_from_baseline[min_idx]
    ) * linewidth_guess_G
    if np.isclose(amp_guess, 0):
        amp_guess = enhancement_span * linewidth_guess_G
    fitdata.set_guess(
        E_0=dict(
            value=baseline_guess,
            min=s.data.real.min() - enhancement_span,
            max=s.data.real.max() + enhancement_span,
        ),
        A=dict(
            value=amp_guess,
            min=-100 * abs(amp_guess),
            max=100 * abs(amp_guess),
        ),
        Bcenter=dict(
            value=center_guess_G, min=B_axis_G.min(), max=B_axis_G.max()
        ),
        sigma=dict(
            value=linewidth_guess_G,
            min=0.01,
            max=max(sweep_width_G, 0.1),
        ),
        lambda_L=dict(
            value=linewidth_guess_G,
            min=0.01,
            max=max(sweep_width_G, 0.1),
        ),
        eta=dict(value=0.5, min=0, max=1),
    )

    # Match the robust path from the earlier script: locate the resonance
    # first, then fit the line shape.
    for name, par in fitdata.guess_parameters.items():
        if name != "Bcenter":
            par.vary = False
    fitdata.fit(use_jacobian=False)
    fitdata.guess_parameters = fitdata.fit_parameters
    for par in fitdata.guess_parameters.values():
        par.vary = True
    fitdata.fit(use_jacobian=False)

    fit = fitdata.eval(800)
    peak_field_G, peak_enhancement = _largest_deviation(
        fit, fitdata.output("E_0")
    )
    return fit, dict(
        model="single_voigt",
        baseline=fitdata.output("E_0"),
        center_field_G=fitdata.output("Bcenter"),
        linewidth_G=_fwhm_pseudo_voigt(
            fitdata.output("sigma"), fitdata.output("lambda_L")
        ),
        peak_field_G=peak_field_G,
        peak_enhancement=peak_enhancement,
        amplitude=fitdata.output("A"),
        gaussian_sigma_G=fitdata.output("sigma"),
        lorentzian_hwhm_G=fitdata.output("lambda_L"),
        lorentzian_fraction=fitdata.output("eta"),
        fit_output=fitdata.output(),
    )


def fit_three_voigt(
    s, center_guesses_G, sigma_guesses_G, lambda_guesses_G, eta_guesses
):
    """
    Fit three absorption pseudo-Voigt lines in gauss.

    Use this for resolved TEMPOL-like spectra.  For arbitrary BSA, peptoid, or
    spin-labeled samples, the single-Voigt model is usually a safer default.
    """
    B, E0 = sp.symbols("B E_0", real=True)
    A0, A1, A2 = sp.symbols("A_0 A_1 A_2", real=True)
    B0, B1, B2 = sp.symbols("B_0 B_1 B_2", real=True)
    sigma0, sigma1, sigma2 = sp.symbols(
        "sigma_0 sigma_1 sigma_2",
        positive=True,
    )
    lambda0, lambda1, lambda2 = sp.symbols(
        "lambda_0 lambda_1 lambda_2",
        positive=True,
    )
    eta0, eta1, eta2 = sp.symbols("eta_0 eta_1 eta_2", real=True)

    fitdata = psd.lmfitdata(s)
    fitdata.functional_form = (
        E0
        + A0 * _pseudo_voigt(B, B0, sigma0, lambda0, eta0)
        + A1 * _pseudo_voigt(B, B1, sigma1, lambda1, eta1)
        + A2 * _pseudo_voigt(B, B2, sigma2, lambda2, eta2)
    )

    B_axis_G = np.asarray(s["B"], dtype=float)
    enhancement_span = max(np.ptp(s.data.real), 0.1)
    amp_guess = s.data.real[np.argmax(abs(s.data.real - 1.0))] - 1.0
    fitdata.set_guess(
        E_0=dict(
            value=1.0,
            min=1.0 - 2 * enhancement_span,
            max=1.0 + 2 * enhancement_span,
        ),
        A_0=dict(
            value=amp_guess,
            min=min(3 * amp_guess, 0),
            max=max(3 * amp_guess, 0),
        ),
        A_1=dict(
            value=amp_guess,
            min=min(3 * amp_guess, 0),
            max=max(3 * amp_guess, 0),
        ),
        A_2=dict(
            value=amp_guess,
            min=min(3 * amp_guess, 0),
            max=max(3 * amp_guess, 0),
        ),
        B_0=dict(
            value=center_guesses_G[0], min=B_axis_G.min(), max=B_axis_G.max()
        ),
        B_1=dict(
            value=center_guesses_G[1], min=B_axis_G.min(), max=B_axis_G.max()
        ),
        B_2=dict(
            value=center_guesses_G[2], min=B_axis_G.min(), max=B_axis_G.max()
        ),
        sigma_0=dict(
            value=sigma_guesses_G[0], min=0.02, max=max(np.ptp(B_axis_G), 1)
        ),
        sigma_1=dict(
            value=sigma_guesses_G[1], min=0.02, max=max(np.ptp(B_axis_G), 1)
        ),
        sigma_2=dict(
            value=sigma_guesses_G[2], min=0.02, max=max(np.ptp(B_axis_G), 1)
        ),
        lambda_0=dict(
            value=lambda_guesses_G[0], min=0.02, max=max(np.ptp(B_axis_G), 1)
        ),
        lambda_1=dict(
            value=lambda_guesses_G[1], min=0.02, max=max(np.ptp(B_axis_G), 1)
        ),
        lambda_2=dict(
            value=lambda_guesses_G[2], min=0.02, max=max(np.ptp(B_axis_G), 1)
        ),
        eta_0=dict(value=eta_guesses[0], min=0, max=1),
        eta_1=dict(value=eta_guesses[1], min=0, max=1),
        eta_2=dict(value=eta_guesses[2], min=0, max=1),
    )

    active_lines = [
        j
        for j, center_G in enumerate(center_guesses_G)
        if B_axis_G.min() - 2 <= center_G <= B_axis_G.max() + 2
    ]
    if not active_lines:
        active_lines = [np.argmin(abs(center_guesses_G - B_axis_G.mean()))]
    for j in set(range(3)) - set(active_lines):
        fitdata.guess_parameters[f"A_{j}"].value = 0

    active_parameters = {"E_0"}
    for j in active_lines:
        active_parameters.update(
            {f"A_{j}", f"B_{j}", f"sigma_{j}", f"lambda_{j}", f"eta_{j}"}
        )

    # First fit only the baseline, amplitudes, and covered centers.  This
    # avoids letting linewidths compensate for poor initial center guesses.
    for name, par in fitdata.guess_parameters.items():
        par.vary = (
            name
            in {
                "E_0",
                "A_0",
                "A_1",
                "A_2",
                "B_0",
                "B_1",
                "B_2",
            }
            & active_parameters
        )
    fitdata.fit(use_jacobian=False)
    fitdata.guess_parameters = fitdata.fit_parameters
    for name, par in fitdata.guess_parameters.items():
        par.vary = name in active_parameters
    fitdata.fit(use_jacobian=False)

    fit = fitdata.eval(800)
    peak_field_G, peak_enhancement = _largest_deviation(
        fit, fitdata.output("E_0")
    )
    line_results = []
    for j in range(3):
        line_results.append(
            dict(
                center_field_G=fitdata.output(f"B_{j}"),
                amplitude=fitdata.output(f"A_{j}"),
                linewidth_G=_fwhm_pseudo_voigt(
                    fitdata.output(f"sigma_{j}"),
                    fitdata.output(f"lambda_{j}"),
                ),
                gaussian_sigma_G=fitdata.output(f"sigma_{j}"),
                lorentzian_hwhm_G=fitdata.output(f"lambda_{j}"),
                lorentzian_fraction=fitdata.output(f"eta_{j}"),
            )
        )
    strongest_line = max(line_results, key=lambda x: abs(x["amplitude"]))
    return fit, dict(
        model="three_voigt",
        baseline=fitdata.output("E_0"),
        center_field_G=strongest_line["center_field_G"],
        linewidth_G=strongest_line["linewidth_G"],
        peak_field_G=peak_field_G,
        peak_enhancement=peak_enhancement,
        lines=line_results,
        active_lines=active_lines,
        fit_output=fitdata.output(),
    )


def _run_fit(s, fit_model, poly_order, conversion):
    "Dispatch to the requested fit model and keep all fits in gauss."
    if fit_model is None:
        return None, dict(
            model=None,
            baseline=1.0,
            center_field_G=None,
            linewidth_G=None,
            peak_field_G=None,
            peak_enhancement=None,
        )
    if fit_model == "single_voigt":
        return fit_single_voigt(s)
    if fit_model == "three_voigt":
        shifted_centers_G = three_voigt_centers_G + (
            conversion["resonance_field_G"] - three_voigt_centers_G[1]
        )
        return fit_three_voigt(
            s,
            shifted_centers_G,
            three_voigt_sigma_G,
            three_voigt_lambda_G,
            three_voigt_eta,
        )
    if fit_model == "polynomial":
        return fit_polynomial(s, order=poly_order)
    raise ValueError(
        "fit_model must be 'single_voigt', 'three_voigt', "
        "'polynomial', or None"
    )


def _add_offset_results(fit_results, conversion):
    "Add EPR-offset versions of the fitted field-domain quantities."
    fit_results = fit_results.copy()
    for field_key, offset_key in (
        ("center_field_G", "center_epr_offset_MHz"),
        ("peak_field_G", "peak_epr_offset_MHz"),
    ):
        if fit_results[field_key] is None:
            fit_results[offset_key] = None
        else:
            fit_results[offset_key] = (
                fit_results[field_key] - conversion["resonance_field_G"]
            ) * conversion["epr_field_conversion_MHz_G"]
    if fit_results["linewidth_G"] is None:
        fit_results["linewidth_MHz"] = None
    else:
        fit_results["linewidth_MHz"] = (
            fit_results["linewidth_G"]
            * conversion["epr_field_conversion_MHz_G"]
        )
    if "lines" in fit_results:
        fit_results["lines"] = [
            dict(
                line,
                center_epr_offset_MHz=(
                    line["center_field_G"] - conversion["resonance_field_G"]
                )
                * conversion["epr_field_conversion_MHz_G"],
                linewidth_MHz=(
                    line["linewidth_G"]
                    * conversion["epr_field_conversion_MHz_G"]
                ),
            )
            for line in fit_results["lines"]
        ]
    return fit_results


s = psd.find_file(
    thisfile,
    exp_type=exp_type,
    expno=nodename,
    lookup=prscr.lookup_table,
)

with psd.figlist_var() as fl:
    fl.basename = thisfile
    acq_params = s.get_prop("acq_params")
    _require_acq_params(acq_params)
    conversion = _calibration_from_acq(acq_params)

    # Use the measured Hall readback as the independent variable.  The NMR
    # carrier/microwave calibration is handled separately below.
    s["indirect"] = _field_readback_G(s)
    s.rename("indirect", "B")
    s.set_units("B", "G")
    s.sort("B")

    # Integrate each phase-cycled echo to one real NMR signal per field point.
    s, _ = prscr.rough_table_of_integrals(
        s,
        fl=fl,
        title=thisfile,
    )
    if s.get_units("B") == "kG":
        s["B"] = np.asarray(s["B"], dtype=float) * 1e3
        s.set_units("B", "G")

    _, edge_signal, left_edge_signal, right_edge_signal = _edge_reference(
        s, n_edge
    )
    s /= edge_signal
    s.name("enhancement")

    fit, fit_results = _run_fit(s, fit_model, poly_order, conversion)
    fit_results = _add_offset_results(fit_results, conversion)
    s_forplot = _to_epr_offset_axis(s, conversion)
    fit_forplot = None if fit is None else _to_epr_offset_axis(fit, conversion)

    fl.next("field-swept DNP enhancement", legend=True)
    ax = plt.gca()
    psd.plot(
        s_forplot,
        "o",
        ax=ax,
        label="experimental enhancement",
    )
    if fit_forplot is not None:
        psd.plot(
            fit_forplot,
            ax=ax,
            label=f"{fit_results['model']} fit",
            alpha=0.8,
        )
    ax.axvline(0, ls=":", color="k", alpha=0.5)
    ax.text(
        0,
        0.96,
        " common resonance field",
        ha="left",
        va="top",
        color="k",
        transform=mpl.transforms.blended_transform_factory(
            ax.transData,
            ax.transAxes,
        ),
    )
    if fit_results["peak_epr_offset_MHz"] is not None:
        ax.axvline(
            fit_results["peak_epr_offset_MHz"],
            ls=":",
            color="C3",
            alpha=0.6,
        )
        ax.text(
            fit_results["peak_epr_offset_MHz"],
            0.80,
            " largest |E-baseline|\n"
            f" {fit_results['peak_epr_offset_MHz']:#0.5g} MHz\n"
            f" E={fit_results['peak_enhancement']:#0.5g}",
            ha="left",
            va="top",
            color="C3",
            transform=mpl.transforms.blended_transform_factory(
                ax.transData,
                ax.transAxes,
            ),
        )
    ax.set_xlabel(
        r"EPR offset, "
        r"$\nu_{\mathrm{EPR}}(B)-\nu_{\mathrm{MW}}$ / MHz"
    )
    ax.set_ylabel(r"$E$")
    ax.grid()

    relative_edge_difference = abs(left_edge_signal - right_edge_signal) / abs(
        edge_signal
    )
    print(
        f"proton carrier frequency: {acq_params['carrierFreq_MHz']:#0.10g} MHz"
    )
    print(
        "microwave resonance frequency: "
        f"{acq_params['uw_dip_center_GHz']:#0.10g} GHz"
    )
    print(
        f"common resonance field: {conversion['resonance_field_G']:#0.10g} G"
    )
    print(
        "nmr_epr_freq_ratio: "
        f"{conversion['nmr_epr_freq_ratio']:#0.10g} MHz/GHz "
        "(numerically ppt)"
    )
    print(
        "epr_field_conversion_MHz_G: "
        f"{conversion['epr_field_conversion_MHz_G']:#0.10g} MHz/G"
    )
    print(f"edge-reference signal: {edge_signal:#0.10g}")
    print(f"left edge reference: {left_edge_signal:#0.10g}")
    print(f"right edge reference: {right_edge_signal:#0.10g}")
    print(f"relative edge mismatch: {relative_edge_difference:#0.5g}")
    print(f"fit model: {fit_results['model']}")
    if fit_results["center_field_G"] is not None:
        print(
            "fit center: "
            f"{fit_results['center_field_G']:#0.10g} G, "
            f"{fit_results['center_epr_offset_MHz']:#0.10g} MHz"
        )
    if fit_results["linewidth_G"] is not None:
        print(
            "fit linewidth: "
            f"{fit_results['linewidth_G']:#0.10g} G, "
            f"{fit_results['linewidth_MHz']:#0.10g} MHz"
        )
    if fit_results["peak_field_G"] is not None:
        print(
            "largest fitted enhancement: "
            f"{fit_results['peak_enhancement']:#0.10g} at "
            f"{fit_results['peak_field_G']:#0.10g} G, "
            f"{fit_results['peak_epr_offset_MHz']:#0.10g} MHz"
        )
    if "lines" in fit_results:
        for j, line in enumerate(fit_results["lines"]):
            print(
                f"line {j}: center "
                f"{line['center_field_G']:#0.10g} G, "
                f"{line['center_epr_offset_MHz']:#0.10g} MHz; "
                f"FWHM {line['linewidth_G']:#0.10g} G"
            )
    if "fit_output" in fit_results:
        print(fit_results["fit_output"])
