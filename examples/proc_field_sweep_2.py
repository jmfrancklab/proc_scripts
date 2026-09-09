"""
hydroxyTEMPO ODNP field sweep: absorptive Voigt spectrum
========================================================

Loads the MW-on and MW-off field-sweep nodes, applies the same phasing and
correlation alignment path used by ``proc_field_sweep.py``, forms
``epsilon = 1 - E`` using the same baseline/sign convention as
``proc_field_sweep.py``, and fits the spectrum on the measured field axis
first with three absorptive Fourier-domain Voigt lines and then with added
dispersive components.

The Voigt fit uses magnetic-field offset in gauss. The displayed x-axis is
converted to the NMR/EPR frequency ratio in MHz/GHz only after fitting.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp
from matplotlib.gridspec import GridSpec
from pyspecProcScripts import (
    L2G,
    correl_align,
    hermitian_function_test,
    select_pathway,
    zeroth_order_ph,
)
from pyspecdata import Q_
from scipy.special import jv

# {{{ changeable parameters

if len(sys.argv) == 4:
    on_file, off_file, exp_type = sys.argv[1:]
elif len(sys.argv) == 6:
    on_file, off_file, exp_type, on_node, off_node = sys.argv[1:]
else:
    on_file, off_file, exp_type = (
        "260818_TMTPDI_field_sweep.h5",
        "260818_TMTPDI_field_sweep.h5",
        "b27/field_dependent",
    )
    on_node, off_node = "field_sweep_1", "field_sweep_2"
signal_pathway = {"ph1": 1}
peak_lower_thresh = 0.1
apodization = True
apod_width_Hz = 1e3
mod_amp = 0.0

# }}}


# SINGLE_USE_EXCEPTION -- control-flow clarity
def phase_node_and_align(filename, nodename, apod=False, fl=None):
    """Load, phase, optionally apodize, and correlation-align one node."""
    s = psd.find_file(
        filename, exp_type=exp_type, expno=nodename, lookup=prscr.lookup_table
    )
    acq_params = s.get_prop("acq_params")
    field_G = s.get_prop("field_readback_G")["NUMPY_DATA"]
    center_field_G = (
        acq_params["carrierFreq_MHz"] / acq_params["gamma_eff_MHz_G"]
    )
    s["indirect"] = (
        (field_G - center_field_G)
        * acq_params["uw_dip_center_GHz"]
        * 1e3
        / center_field_G
    )
    s.rename("indirect", "nu_offset").set_units("nu_offset", "MHz")
    s.sort("nu_offset")
    s.set_prop("coherence_pathway", signal_pathway)

    if apod:
        fig = fl.next(f"{nodename} apodized signal")
        gs = GridSpec(1, 3, figure=fig, left=0.05, right=0.95)
        psd.DCCT(s, fig, title="Raw Data", bbox=gs[0, 0])
        s.ift("t2")
        s *= L2G(apod_width_Hz, criterion="energy")(s.fromaxis("t2"))
        s.ft("t2")
        psd.DCCT(
            s,
            fig,
            title="Equal Energy Apodization",
            bbox=gs[0, 1],
        )

    signal = select_pathway(s, signal_pathway)
    nu_axis = signal.getaxis("nu_offset")
    frq_center, frq_half = prscr.find_peakrange(
        signal["nu_offset" : nu_axis[abs(nu_axis).argmin()]],
        peak_lower_thresh=peak_lower_thresh,
    )
    frq_half = abs(frq_half)

    s.ift("t2")
    s /= zeroth_order_ph(select_pathway(s, signal_pathway))
    s["t2"] -= s.getaxis("t2")[0]
    best_shift = hermitian_function_test(
        select_pathway(s.C.mean("nu_offset"), signal_pathway)
    )
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    s.ft("t2")

    def frq_mask(this_s):
        repeat_dim = (
            "repeats" if "repeats" in this_s.dimlabels else "nu_offset"
        )
        frq_center_local = (
            select_pathway(this_s, signal_pathway)
            .mean(repeat_dim)
            .argmax("t2")
        )
        return this_s * np.exp(
            -((this_s.fromaxis("t2") - frq_center_local) ** 2)
            / (4 * frq_half**2)
        )

    def coherence_unmask_fn(coh_array):
        thisslice = coh_array
        for j, (k, v) in enumerate(signal_pathway.items()):
            if j == len(signal_pathway) - 1:
                thisslice[k, v] = 1
            else:
                thisslice = thisslice[k, v]
        return coh_array

    repeat_sign = select_pathway(s, signal_pathway).C.real.sum("t2")
    repeat_sign = repeat_sign.run(np.sign)
    opt_shift = correl_align(
        s * repeat_sign,
        frq_mask_fn=frq_mask,
        coherence_unmask_fn=coherence_unmask_fn,
        repeat_dims="nu_offset",
        max_shift=frq_half,
    )
    s.ift("t2").ift(list(signal_pathway))
    s *= np.exp(-1j * 2 * np.pi * opt_shift * s.fromaxis("t2"))
    s.ft("t2").ft(list(signal_pathway.keys()))

    if apod:
        psd.DCCT(
            s,
            fig,
            title="After apodization\nand alignment",
            bbox=gs[0, 2],
        )

    return s, acq_params, center_field_G, frq_center, frq_half


# SINGLE_USE_EXCEPTION -- control-flow clarity
def fit_in_stages(fitdata, stages):
    """Run staged lmfit optimization using current guesses as warm starts."""
    for stage_number, (label, vary_prefixes) in enumerate(stages):
        if stage_number > 0:
            fitdata.guess_parameters = fitdata.fit_parameters
        for name, parameter in fitdata.guess_parameters.items():
            parameter.vary = any(
                name.startswith(prefix) for prefix in vary_prefixes
            )
        print(f"about to run {label}")
        fitdata.fit(use_jacobian=False)


with psd.figlist_var() as fl:
    fl.basename = on_file

    (
        off,
        _,
        _,
        off_frq_center,
        off_frq_half,
    ) = phase_node_and_align(off_file, off_node, apod=apodization, fl=fl)
    (
        on,
        acq_params,
        center_field_G,
        on_frq_center,
        on_frq_half,
    ) = phase_node_and_align(
        on_file,
        on_node,
        apod=apodization,
        fl=fl,
    )

    on_band = select_pathway(on, signal_pathway)[
        "t2" : (on_frq_center - on_frq_half, on_frq_center + on_frq_half)
    ].integrate("t2")
    off_band = select_pathway(off, signal_pathway)[
        "t2" : (off_frq_center - off_frq_half, off_frq_center + off_frq_half)
    ].integrate("t2")

    nu_on = np.asarray(on_band["nu_offset"], dtype=float)
    epr_shift_per_G = acq_params["uw_dip_center_GHz"] * 1e3 / center_field_G
    # Thermal signal should be common across the sweep, so use one positive
    # flat baseline averaged over all MW-off field points.
    thermal_level = abs(np.mean(off_band.data.real))
    field_offset_on = nu_on / epr_shift_per_G
    field_spacing_guess = 44.5 / epr_shift_per_G
    expected_centers = field_spacing_guess * np.r_[-1, 0, 1]
    peak_mask = np.zeros_like(field_offset_on, dtype=bool)
    for center_guess in expected_centers:
        peak_mask |= (
            abs(field_offset_on - center_guess) < 0.45 * field_spacing_guess
        )
    if not peak_mask.any():
        peak_mask[:] = True

    enhancement_for_sign = on_band.data.real / thermal_level
    peak_enhancement = enhancement_for_sign[peak_mask]
    strongest_peak = peak_enhancement[np.argmax(abs(peak_enhancement))]
    if strongest_peak > 0:
        on_band *= -1
    epsilon = on_band.C
    epsilon.data = 1 - on_band.data.real / thermal_level
    epsilon.set_error(None)
    epsilon.name("epsilon")

    epsilon.rename("nu_offset", "B")
    epsilon.setaxis("B", field_offset_on).set_units("B", "G")

    field_fit = np.linspace(
        field_offset_on.min(), field_offset_on.max(), len(field_offset_on)
    )
    fit_epsilon = psd.nddata(
        np.interp(field_fit, field_offset_on, epsilon.data.real),
        ["B"],
    )
    fit_epsilon.labels("B", field_fit)
    fit_epsilon.name("epsilon")

    edge_count = max(3, len(field_offset_on) // 10)
    edge_baseline = np.r_[
        fit_epsilon.data.real[:edge_count],
        fit_epsilon.data.real[-edge_count:],
    ].mean()
    peak_signal = fit_epsilon.C
    peak_signal -= edge_baseline
    signal_scale = float(peak_signal.data.real.max())
    if signal_scale <= 0:
        raise ValueError("Cannot find positive peaks in epsilon")

    peak_centers = []
    peak_heights = []
    for center_guess in expected_centers:
        local_signal = peak_signal[
            "B" : (
                center_guess - 0.45 * field_spacing_guess,
                center_guess + 0.45 * field_spacing_guess,
            )
        ]
        if local_signal.data.size == 0:
            peak_centers.append(center_guess)
            peak_heights.append(peak_signal["B":center_guess].item().real)
        else:
            peak_centers.append(center_guess)
            peak_heights.append(local_signal.data.real.max())
    linewidth_guess = field_spacing_guess / 2

    B = sp.symbols("B", real=True)
    A, E_0, lambda_L, Bcenter, sigma = sp.symbols(
        "A E_0 lambda_L Bcenter sigma", real=True
    )
    voigt_line = (
        A
        * sp.exp(1j * 2 * np.pi * Bcenter * B)
        * sp.exp(
            -lambda_L * sp.pi * abs(B) - sp.pi**2 * abs(B) ** 2 * sigma**2
        )
    )
    A_symbols = sp.symbols("A0:3", real=True)
    A_disp_symbols = sp.symbols("A_disp0:3", real=True)
    Bcenter_symbols = sp.symbols("Bcenter0:3", real=True)
    FWHM_symbols = sp.symbols("FWHM0:3", real=True)
    L_vs_G_frac_symbols = sp.symbols("L_vs_G_frac0:3", real=True)
    voigt_fwhm_coeff = 0.5346
    voigt_fwhm_coeff_symb = sp.Float(voigt_fwhm_coeff)
    voigt_fwhm_remainder_symb = (sp.Integer(1) - voigt_fwhm_coeff_symb) ** 2

    def build_model(amplitudes):
        thefunction = 0
        for j, amplitude in enumerate(amplitudes):
            lorentzian_FWHM = FWHM_symbols[j] * L_vs_G_frac_symbols[j]
            gaussian_FWHM = sp.sqrt(
                (FWHM_symbols[j] - voigt_fwhm_coeff_symb * lorentzian_FWHM)
                ** 2
                - voigt_fwhm_remainder_symb * lorentzian_FWHM**2
            )
            thefunction += voigt_line.subs(
                {A: amplitude, Bcenter: Bcenter_symbols[j]}
            ).subs(
                {
                    lambda_L: lorentzian_FWHM,
                    sigma: gaussian_FWHM / (2 * sp.sqrt(sp.log(2))),
                }
            )
        return thefunction

    fit_input = fit_epsilon.C
    fit_input.ift("B", shift=True)
    fitdata = psd.lmfitdata(fit_input)

    @fitdata.define_data_transform
    def my_data_transform(d_local):
        d_local.ft("B")
        return d_local.real

    @fitdata.define_residual_transform
    def my_residual_transform(d_local):
        h_m = Q_(mod_amp, "G").to("G").magnitude
        d_local *= d_local.fromaxis(
            "B",
            lambda axis: (
                jv(0, abs(h_m * axis * np.pi)) - jv(2, abs(h_m * axis * np.pi))
            ),
        )
        d_local.ft("B")
        return d_local.real

    # A constant spectral baseline is a zero-time impulse in the FT domain.
    baseline_impulse = fit_input.C
    baseline_impulse.data[:] = 0
    baseline_impulse.data[abs(baseline_impulse.getaxis("B")).argmin()] = 1
    baseline_transform_scale = my_residual_transform(
        baseline_impulse
    ).data.real.mean()
    fitdata.functional_form = build_model(
        A_symbols
    ) + E_0 / baseline_transform_scale * sp.DiracDelta(B)

    fit_guesses = {
        "E_0": {
            "value": edge_baseline,
            "min": edge_baseline - signal_scale,
            "max": edge_baseline + signal_scale,
        }
    }
    for j, center in enumerate(peak_centers):
        area_guess = max(peak_heights[j], 0.1 * signal_scale)
        area_guess *= linewidth_guess
        fit_guesses[f"A{j}"] = {
            "value": area_guess,
            "min": 0,
            "max": 20 * area_guess,
        }
        fit_guesses[f"Bcenter{j}"] = {
            "value": center,
            "min": expected_centers[j] - 0.35 * field_spacing_guess,
            "max": expected_centers[j] + 0.35 * field_spacing_guess,
        }
        fit_guesses[f"FWHM{j}"] = {
            "value": linewidth_guess,
            "min": 0.1 * linewidth_guess,
            "max": 2.6 * linewidth_guess,
        }
        fit_guesses[f"L_vs_G_frac{j}"] = {
            "value": 0.5,
            "min": 0.01,
            "max": 0.99,
        }
    fitdata.set_guess(fit_guesses)
    fit_in_stages(
        fitdata,
        [
            ("vary only centers", ("Bcenter",)),
            ("vary centers and amplitudes", ("Bcenter", "A", "E_0")),
            (
                "vary centers, amplitudes, and FWHM",
                ("Bcenter", "A", "E_0", "FWHM"),
            ),
            ("final absorptive Voigt stage", tuple(fitdata.guess_parameters)),
        ],
    )

    display_range = (field_offset_on.min(), field_offset_on.max())
    absorptive_curve = fitdata.eval(500)["B":display_range]
    absorptive_parameters = fitdata.fit_parameters.copy()
    absorptive_output = fitdata.output()
    absorptive_redchi = fitdata.fit_output.redchi

    fitdata.functional_form = (
        build_model(A_symbols)
        - sp.I * sp.sign(B) * build_model(A_disp_symbols)
        + E_0 / baseline_transform_scale * sp.DiracDelta(B)
    )
    extended_guesses = {
        name: {
            "value": np.clip(
                parameter.value,
                parameter.min + 0.02 * (parameter.max - parameter.min),
                parameter.max - 0.02 * (parameter.max - parameter.min),
            ),
            "min": parameter.min,
            "max": parameter.max,
        }
        for name, parameter in absorptive_parameters.items()
        if name in fitdata.guess_parameters
    }
    for j in range(3):
        absorptive_area = absorptive_parameters[f"A{j}"].value
        dispersive_guess = 0
        if j == 0:
            dispersive_guess = -absorptive_area
        elif j == 2:
            dispersive_guess = absorptive_area
        extended_guesses[f"A_disp{j}"] = {
            "value": dispersive_guess,
            "min": -10 * absorptive_area,
            "max": 10 * absorptive_area,
        }
    fitdata.set_guess(extended_guesses)
    fit_in_stages(
        fitdata,
        [
            (
                "center fields and all amplitudes",
                ("A", "Bcenter", "E_0"),
            ),
            (
                "full absorptive/dispersive stage",
                tuple(fitdata.guess_parameters),
            ),
        ],
    )
    fitted_curve = fitdata.eval(500)["B":display_range]
    fit_output = fitdata.output()
    component_curves = []
    for j in range(3):
        component_guesses = fit_output.copy()
        for k in range(3):
            if k != j:
                component_guesses[f"A{k}"] = 0
                component_guesses[f"A_disp{k}"] = 0
        component_guesses["E_0"] = 0
        fitdata.set_guess(component_guesses)
        component_curves.append(
            fitdata.set_to_guess().eval(500)["B":display_range]
        )
    fitdata.set_guess(fit_output)
    for curve in [
        epsilon,
        absorptive_curve,
        fitted_curve,
        *component_curves,
    ]:
        curve.setaxis(
            "B",
            lambda field_offset: (
                (center_field_G + field_offset)
                * acq_params["gamma_eff_MHz_G"]
                / acq_params["uw_dip_center_GHz"]
            ),
        )
        curve.rename("B", "MHz_per_GHz")
        curve.set_units("MHz_per_GHz", "MHz/GHz")

    fig, ax = plt.subplots(figsize=(4, 4))
    fl.next("hydroxyTEMPO absorptive and dispersive ODNP field sweep", fig=fig)
    psd.plot(
        epsilon,
        "o",
        ax=ax,
        color="k",
        markersize=3,
        label=r"$\epsilon = 1 - E$",
    )
    psd.plot(
        absorptive_curve,
        ax=ax,
        color="C1",
        linestyle="--",
        linewidth=1.5,
        label="absorptive Voigt fit",
    )
    psd.plot(
        fitted_curve,
        ax=ax,
        color="r",
        linewidth=2,
        label="absorptive + dispersive Voigt fit",
    )
    for j, component in enumerate(component_curves):
        psd.plot(
            component,
            ax=ax,
            linestyle="--",
            linewidth=1.5,
            label=f"component {j + 1}",
        )
    line_centers = tuple(
        (center_field_G + fit_output[f"Bcenter{j}"])
        * acq_params["gamma_eff_MHz_G"]
        / acq_params["uw_dip_center_GHz"]
        for j in range(3)
    )
    ratio_center = line_centers[1]
    ratio_spacing = 0.5 * (line_centers[2] - line_centers[0])
    for this_center in line_centers:
        ax.axvline(this_center, ls=":", color="grey", alpha=0.5)
    ax.axhline(0, lw=0.6, color="k")
    ax.set_xlim(
        epsilon.getaxis("MHz_per_GHz").min(),
        epsilon.getaxis("MHz_per_GHz").max(),
    )
    ax.set_xlabel("MHz/GHz")
    ax.set_ylabel(r"$\epsilon = 1 - E$")
    ax.grid(alpha=0.3)
    fig.tight_layout()

    print(f"center field: {center_field_G:#0.8g} G")
    print(f"peak epsilon: {epsilon.data.max():#0.4g}")
    print(f"absorptive fit reduced chi-square: {absorptive_redchi:#0.5g}")
    print(
        "absorptive + dispersive fit reduced chi-square: "
        f"{fitdata.fit_output.redchi:#0.5g}"
    )
    print(
        "absorptive line centers (MHz/GHz): "
        + ", ".join(
            "{:#0.8g}".format(
                (center_field_G + absorptive_output[f"Bcenter{j}"])
                * acq_params["gamma_eff_MHz_G"]
                / acq_params["uw_dip_center_GHz"]
            )
            for j in range(3)
        )
    )
    print(f"NMR/EPR ratio from middle peak: {ratio_center:#0.8g} MHz/GHz")
    print(f"fit ratio spacing: {ratio_spacing:#0.4g} MHz/GHz")
    print(
        "line centers (MHz/GHz): "
        + ", ".join(f"{center:#0.8g}" for center in line_centers)
    )
