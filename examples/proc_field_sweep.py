"""
hydroxyTEMPO ODNP field sweep: DNP spectrum (epsilon = 1 - E)
==============================================================

``py proc_field_sweep.py ON_FILE OFF_FILE EXP_TYPE [ON_NODE OFF_NODE]``

Loads the MW-on and MW-off field-sweep files/nodes with
``find_file(..., lookup=lookup_table)`` (dispatching on the saved
``field_sweep_v5`` postproc, so each node returns conjugated,
``tau``-shifted, FT'd along ``t2`` and along the phase cycle, with
``coherence_pathway`` / ``acq_params`` / ``field_readback_G`` as
properties), phases each node, integrates the on-resonance NMR peak for
every field, and plots ``epsilon = 1 - E`` (``E = S_on / S_off``) with a
three-line Voigt fit via ``lmfitdata``, alongside the analytic derivative
of the fit.

Phasing (before any integration or plotting)
---------------------------------------------
Per repo convention (see ``examples/Hermitian_Phasing.py``), each node gets:

Opt. **Equal-energy apodization**
1. **zeroth-order phasing** (``zeroth_order_ph``) -- a single global phase
   for the whole node, found from the moment of inertia of the
   coherence-selected signal in the complex plane.
2. **Hermitian (first-order) phasing** (``hermitian_function_test``) -- a
   single global timing correction that centers the echo, found from the
   nu_offset-averaged, coherence-selected FID/echo.
3. **Correlation alignment**

Sign convention
----------------
``E = S_on / S_off`` is negative under the DNP lines (1H Overhauser
enhancement inverts the water-proton signal).  We plot the more
conventional ``epsilon = 1 - E`` instead, which is positive under the DNP
lines and zero off-resonance.
"""

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
from scipy.special import wofz

# {{{ changeable parameters

on_file, off_file, exp_type = (
    "260720_hydroxytempo_field_sweep.h5",
    "260720_hydroxytempo_field_sweep.h5",
    "b27/field_dependent",
)
on_node, off_node = "field_sweep_1", "field_sweep_2"  # MW on / MW off
signal_pathway = {"ph1": 1}
peak_lower_thresh = 0.1
apodization = True  # Apply to both or none. Recommended.
apod_width_Hz = 1e3  # After the optimal value increasing shouldn't change
# the enhancement. Low values give weird results.
# }}}}


def phase_node_and_align(filename, nodename, apod=False, fl=None):
    """Load one node and apply zeroth-order phasing and Hermitian timing
    correction.  Optionally apply equal-energy apodization and correlation
    alignment.  Returns the phased nddata, acq_params, the center field, and
    the zeroth-order phase actually used (so the caller can share it with the
    other node)."""
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
        # {{{ Equal-energy apodization
        fig = fl.next(f"{nodename} apodized signal")
        fl.title = f"{nodename} apodized signal"
        gs = GridSpec(1, 3, figure=fig, left=0.05, right=0.95)
        psd.DCCT(
            s,
            fig,
            title="Raw Data",
            bbox=gs[0, 0],
        )
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
    zeroth_phase = zeroth_order_ph(select_pathway(s, signal_pathway))
    s /= zeroth_phase
    s["t2"] -= s.getaxis("t2")[0]  # needed for hermitian_function_test
    best_shift = hermitian_function_test(
        select_pathway(s.C.mean("nu_offset"), signal_pathway)
    )
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    s.ft("t2")

    def frq_mask(this_s):
        repeat_dim = (
            "repeats" if "repeats" in this_s.dimlabels else "nu_offset"
        )
        nu_center = (
            select_pathway(this_s, signal_pathway)
            .mean(repeat_dim)
            .argmax("t2")
        )
        return this_s * np.exp(
            -((this_s.fromaxis("t2") - nu_center) ** 2) / (4 * frq_half**2)
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


with psd.figlist_var() as fl:
    fl.basename = on_file

    (
        off,
        acq_params,
        center_field_G,
        off_frq_center,
        off_frq_half,
    ) = phase_node_and_align(off_file, off_node, apod=apodization, fl=fl)
    on, _, _, on_frq_center, on_frq_half = phase_node_and_align(
        on_file, on_node, apod=apodization, fl=fl
    )

    on_band = select_pathway(on, signal_pathway)[
        "t2" : (on_frq_center - on_frq_half, on_frq_center + on_frq_half)
    ].integrate("t2")
    off_band = select_pathway(off, signal_pathway)[
        "t2" : (off_frq_center - off_frq_half, off_frq_center + off_frq_half)
    ].integrate("t2")

    # {{{ epsilon = 1 - E, E = S_on / S_off, off node interpolated onto the
    #     on-node field axis.
    nu_on = np.asarray(on_band["nu_offset"], dtype=float)
    nu_off = np.asarray(off_band["nu_offset"], dtype=float)
    off_order = np.argsort(nu_off)
    off_on_axis = np.interp(
        nu_on, nu_off[off_order], off_band.data.real[off_order]
    )
    thermal_level = np.median(off_on_axis)
    # DNP data is negative wrt thermals
    if np.median(on_band.data.real) * thermal_level > 0:
        on_band *= -1
    epsilon = on_band.C
    epsilon.data = 1 - on_band.data.real / thermal_level
    epsilon.set_error(None)  # None for now. In the future we can add errors
    # from the goodness of the fit.

    epsilon.name("epsilon")
    # }}}

    # {{{ Three-line Voigt fit.
    nu_offset, E_0 = sp.symbols("nu_offset E_0", real=True)
    A_symbols = sp.symbols("A0:3", real=True)
    A_disp_symbols = sp.symbols("A_disp0:3", real=True)
    Bcenter_symbols = sp.symbols("Bcenter0:3", real=True)
    FWHM_symbols = sp.symbols("FWHM0:3", real=True)
    L_vs_G_frac_symbols = sp.symbols("L_vs_G_frac0:3", real=True)
    wofz_symbol = sp.Function("wofz")
    voigt_fwhm_coeff = 0.5346
    voigt_fwhm_remainder = (1.0 - voigt_fwhm_coeff) ** 2

    def build_model(absorptive_amplitudes, dispersive_amplitudes=None):
        model = E_0
        for j, absorptive_amplitude in enumerate(absorptive_amplitudes):
            lorentzian_FWHM = FWHM_symbols[j] * L_vs_G_frac_symbols[j]
            gaussian_FWHM = sp.sqrt(
                (
                    FWHM_symbols[j]
                    - sp.Float(voigt_fwhm_coeff) * lorentzian_FWHM
                )
                ** 2
                - sp.Float(voigt_fwhm_remainder) * lorentzian_FWHM**2
            )
            gaussian_sigma = gaussian_FWHM / (2 * sp.sqrt(2 * sp.log(2)))
            z = (
                nu_offset - Bcenter_symbols[j] + sp.I * lorentzian_FWHM / 2
            ) / (gaussian_sigma * sp.sqrt(2))
            complex_voigt = wofz_symbol(z) / (
                gaussian_sigma * sp.sqrt(2 * sp.pi)
            )
            model += absorptive_amplitude * sp.re(complex_voigt)
            if dispersive_amplitudes is not None:
                model += dispersive_amplitudes[j] * sp.im(complex_voigt)
        return model

    def fit_in_stages(stages):
        for stage_number, (label, vary_prefixes) in enumerate(stages):
            if stage_number > 0:
                fitdata.guess_parameters = fitdata.fit_parameters
            for name, parameter in fitdata.guess_parameters.items():
                parameter.vary = any(
                    name.startswith(prefix) for prefix in vary_prefixes
                )
            print(f"about to run {label}")
            fitdata.fit(use_jacobian=False)

    edge_count = max(3, len(nu_on) // 10)
    edge_baseline = np.r_[
        epsilon.data.real[:edge_count],
        epsilon.data.real[-edge_count:],
    ].mean()
    peak_signal = epsilon.C
    peak_signal -= edge_baseline
    signal_scale = float(peak_signal.data.real.max())
    if signal_scale <= 0:
        raise ValueError("Cannot find positive peaks in the enhancement")
    for threshold in np.linspace(0.80, 0.10, 15):
        cutoff = threshold * signal_scale
        peak_ranges = peak_signal.contiguous(
            lambda values, cutoff=cutoff: values.real > cutoff
        )
        peak_chunks = sorted(
            [
                tuple(peak_ranges[j, :])
                for j in range(peak_ranges.shape[0])
                if abs(peak_ranges[j, 1] - peak_ranges[j, 0]) > 0
            ],
            key=lambda chunk: 0.5 * (chunk[0] + chunk[1]),
        )
        if len(peak_chunks) == 3:
            break
    else:
        raise ValueError("Could not find three enhancement peaks")

    peak_centers = [
        peak_signal["nu_offset":chunk].argmax("nu_offset").item()
        for chunk in peak_chunks
    ]
    hyperfine_guess = float(np.median(np.diff(peak_centers)))
    linewidth_guess = hyperfine_guess / 2
    peak_guesses = []
    for center in peak_centers:
        peak_height = peak_signal["nu_offset":center].item().real
        area_guess = max(peak_height, 0.1 * signal_scale) * linewidth_guess
        peak_guesses.append(
            {
                "A": {
                    "value": area_guess,
                    "min": 0,
                    "max": 10 * area_guess,
                },
                "FWHM": {
                    "value": linewidth_guess,
                    "min": 0.25 * linewidth_guess,
                    "max": 4 * linewidth_guess,
                },
                "L_vs_G_frac": {
                    "value": 0.5,
                    "min": 0.01,
                    "max": 0.99,
                },
                "Bcenter": {
                    "value": center,
                    "min": center - linewidth_guess,
                    "max": center + linewidth_guess,
                },
            }
        )

    fitdata = psd.lmfitdata(epsilon)
    fitdata.functional_form = build_model(A_symbols)
    fit_guesses = {
        "E_0": {
            "value": edge_baseline,
            "min": edge_baseline - signal_scale,
            "max": edge_baseline + signal_scale,
        }
    }
    for j, line_guess in enumerate(peak_guesses):
        for name in ("A", "FWHM", "L_vs_G_frac", "Bcenter"):
            fit_guesses[f"{name}{j}"] = line_guess[name]
    fitdata.set_guess(fit_guesses)
    fit_in_stages(
        [
            ("vary only center fields", ("Bcenter",)),
            (
                "vary center fields and amplitudes",
                ("Bcenter", "A", "E_0"),
            ),
            (
                "vary center fields, amplitudes, and FWHM",
                ("Bcenter", "A", "E_0", "FWHM"),
            ),
            (
                "three-Voigt final stage",
                tuple(fitdata.guess_parameters.keys()),
            ),
        ]
    )
    three_voigt_parameters = fitdata.fit_parameters.copy()
    three_voigt_fit = fitdata.eval(500)

    # Add the dispersive components only after the absorptive fit converges.
    fitdata.functional_form = build_model(A_symbols, A_disp_symbols)
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
        for name, parameter in three_voigt_parameters.items()
        if name in fitdata.guess_parameters
    }
    for j in range(3):
        absorptive_area = three_voigt_parameters[f"A{j}"].value
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
        [
            (
                "center fields and all amplitudes",
                ("A", "Bcenter", "E_0"),
            ),
            (
                "full absorptive/dispersive stage",
                tuple(fitdata.guess_parameters.keys()),
            ),
        ]
    )
    fitted_curve = fitdata.eval(500)
    fit_output = fitdata.output()
    # }}}

    # {{{ Analytic derivative of the fitted absorption/dispersive curve.
    # For w(z), dw/dz = -2*z*w(z) + 2i/sqrt(pi).
    nu_fine = np.linspace(nu_on.min(), nu_on.max(), 500)
    derivative_values = np.zeros_like(nu_fine)
    for j in range(3):
        FWHM_out = fit_output[f"FWHM{j}"]
        L_vs_G_frac_out = fit_output[f"L_vs_G_frac{j}"]
        lorentzian_FWHM = FWHM_out * L_vs_G_frac_out
        gaussian_FWHM = np.sqrt(
            max(
                (FWHM_out - voigt_fwhm_coeff * lorentzian_FWHM) ** 2
                - voigt_fwhm_remainder * lorentzian_FWHM**2,
                0,
            )
        )
        gaussian_sigma = gaussian_FWHM / (2 * np.sqrt(2 * np.log(2)))
        z = (
            nu_fine - fit_output[f"Bcenter{j}"] + 1j * lorentzian_FWHM / 2
        ) / (gaussian_sigma * np.sqrt(2))
        profile_derivative = (-2 * z * wofz(z) + 2j / np.sqrt(np.pi)) / (
            2 * gaussian_sigma**2 * np.sqrt(np.pi)
        )
        derivative_values += (
            fit_output[f"A{j}"] * profile_derivative.real
            + fit_output[f"A_disp{j}"] * profile_derivative.imag
        )
    derivative = psd.nddata(derivative_values, ["nu_offset"])
    derivative.setaxis("nu_offset", nu_fine).set_units("nu_offset", "MHz")
    derivative.name("d(epsilon)/d(nu_offset)")
    # }}}

    fig, (ax_eps, ax_deriv) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))
    fl.next("hydroxyTEMPO ODNP field sweep", fig=fig)
    psd.plot(
        epsilon,
        "o",
        ax=ax_eps,
        color="k",
        markersize=3,
        label="DNP curve",
    )
    psd.plot(
        three_voigt_fit,
        ax=ax_eps,
        color="C1",
        linestyle="--",
        alpha=0.6,
        label="3-Voigt fit",
    )
    psd.plot(
        fitted_curve,
        ax=ax_eps,
        color="r",
        alpha=0.8,
        label="3-Voigt + dispersive fit",
    )
    line_centers = tuple(fit_output[f"Bcenter{j}"] for j in range(3))
    nu_c_fit = line_centers[1]
    A_hf_fit = 0.5 * (line_centers[2] - line_centers[0])
    epr_shift_per_G = acq_params["uw_dip_center_GHz"] * 1e3 / center_field_G
    middle_peak_nmr_MHz = (
        center_field_G + nu_c_fit / epr_shift_per_G
    ) * acq_params["gamma_eff_MHz_G"]
    for this_center in line_centers:
        ax_eps.axvline(this_center, ls=":", color="grey", alpha=0.5)
    ax_eps.axhline(0, lw=0.6, color="k")
    ax_eps.set_ylabel(r"$\epsilon  = 1 - E$")
    ax_eps.legend()
    ax_eps.grid(alpha=0.3)

    psd.plot(
        derivative,
        ax=ax_deriv,
        color="b",
        label=r"$d\epsilon/d\nu$",
    )
    for this_center in line_centers:
        ax_deriv.axvline(this_center, ls=":", color="grey", alpha=0.5)
    ax_deriv.axhline(0, lw=0.6, color="k")
    ax_deriv.set_xlabel(r"EPR offset $\nu-\nu_\mathrm{center}$ / MHz")
    ax_deriv.set_ylabel("First Derivative Spectrum")
    ax_deriv.legend()
    ax_deriv.grid(alpha=0.3)
    fig.tight_layout()

    print(f"center field: {center_field_G:#0.8g} G")
    print(f"peak epsilon: {epsilon.data.max():#0.4g}")
    print(
        f"14N hyperfine spacing: {A_hf_fit:#0.4g} MHz; centers (MHz): "
        + ", ".join(f"{center:#0.4g}" for center in line_centers)
    )
    for j in range(3):
        print(
            f"line {j + 1}: FWHM {fit_output[f'FWHM{j}']:#0.4g} MHz, "
            f"L/G {fit_output[f'L_vs_G_frac{j}']:#0.3g}, "
            f"A_disp/A "
            f"{fit_output[f'A_disp{j}'] / fit_output[f'A{j}']:#0.3g}"
        )
    print(
        "NMR/EPR ratio: "
        f"{middle_peak_nmr_MHz / acq_params['uw_dip_center_GHz']:#0.8g}"
        " MHz/GHz "
    )
