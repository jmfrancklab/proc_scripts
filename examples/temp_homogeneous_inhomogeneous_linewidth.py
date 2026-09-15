r"""
Homogeneous and Inhomogeneous Linewidth Determination
=====================================================

For the same experimental FIR and ODNP-enhancement acquisitions, show that the
homogeneous and inhomogeneous linewidths are independent measurements.  The
homogeneous Lorentzian linewidth is fitted from the centered, causal FID-side
decay.  The inhomogeneous linewidth is the separation between the outer edges
of the valid 3%-of-maximum region in the observed ensemble spectrum.

Equal-energy L2G apodization is shown only as the coarse low-SNR peak locator.
It does not set the final inhomogeneous bounds, which come from the original
spectrum with only the small acquisition-length conditioning applied.

This is a temporary diagnostic example and is not intended for committing.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

import pyspecProcScripts as pypcs
from temp_bounds_dcct_data import (
    EXPERIMENT_TYPES,
    determine_preliminary_bounds,
    load_dcct_dataset,
    show_raw_data,
    sign_corrected,
)

plt.rcParams["image.aspect"] = "auto"
with psd.figlist_var() as fl:
    for experiment_type in EXPERIMENT_TYPES:
        raw_data, _, configuration = load_dcct_dataset(experiment_type)
        show_raw_data(fl, raw_data, configuration)

        # Use the narrow preliminary range only to remove known FIR or
        # enhancement sign changes.  The final inhomogeneous bounds below are
        # independently redetermined at the 3% threshold.
        linewidth_data = sign_corrected(
            raw_data,
            configuration,
            determine_preliminary_bounds(raw_data, configuration)[
                "inh_bounds"
            ],
        )
        echo_center = pypcs.find_exponential_echo_center(
            linewidth_data,
            decay_rate=250,
        )
        linewidth_data.set_prop("echo_center", echo_center)
        fid_side = pypcs.fid_side_from_echo(
            linewidth_data,
            echo_center,
        )
        homogeneous_fit_name = f"{experiment_type}: homogeneous linewidth fit"
        homogeneous_linewidth, used_lsq_fallback = pypcs.fit_envelope(
            pypcs.select_pathway(
                fid_side,
                configuration["signal_pathway"],
            ),
            plot_name=homogeneous_fit_name,
            mult_two=True,
            return_fallback=True,
            fl=fl,
        )
        fl.next(homogeneous_fit_name)
        # fit_envelope includes the dimensionless apodization window on its
        # diagnostic; its default data-name label otherwise looks like data.
        plt.gca().lines[-1].set_label("equal-energy L2G window")
        plt.gca().set_title(
            f"{experiment_type}: homogeneous linewidth "
            f"{homogeneous_linewidth:.1f} Hz"
            + (" (provisional LSQ)" if used_lsq_fallback else "")
        )
        # The acquisition is much longer than the fitted decay.  Restrict the
        # diagnostic to ten decay constants so the signal-bearing interval is
        # visible instead of compressing it into the first few pixels.
        plt.gca().set_xlim(0, 10 / (np.pi * homogeneous_linewidth) * 1e3)

        pypcs.det_inh_bounds(
            linewidth_data,
            0.10,
            peak_lowest_thresh=0.03,
            signal_pathway=configuration["signal_pathway"],
            apodization_linewidth=homogeneous_linewidth,
        )
        inhomogeneous_bounds = linewidth_data.get_prop("inh_bounds").copy()
        inhomogeneous_linewidth = np.diff(inhomogeneous_bounds).item()

        # Reproduce the detector spectra for a transparent visualization of
        # the thresholds.  Only this copy receives the weak 1/(5 aq)
        # conditioning; neither the linewidth fit nor the raw data is changed.
        detector_fid = fid_side.C
        acquisition_time = 1 / abs(detector_fid.get_ft_prop("t2", "df"))
        detector_fid *= np.exp(
            -abs(detector_fid.fromaxis("t2")) / (5 * acquisition_time)
        )
        detector_fid = pypcs.select_pathway(
            detector_fid,
            configuration["signal_pathway"],
        )
        l2g_detector_fid = detector_fid.C
        l2g_detector_fid *= pypcs.L2G(
            homogeneous_linewidth,
            criterion="energy",
        )(l2g_detector_fid.fromaxis("t2"))

        detector_spectra = []
        for detector_input in [detector_fid, l2g_detector_fid]:
            detector_spectrum = detector_input.C.ft("t2")
            detector_spectrum.mean_all_but(["t2"]).run(abs)
            raw_detector_spectrum = detector_spectrum.C
            detector_spectrum.convolve(
                "t2",
                50.0,
                enforce_causality=False,
            )
            spectral_width = 1 / detector_spectrum.get_ft_prop("t2", "dt")
            detector_spectrum -= (
                detector_spectrum[
                    "t2" : tuple(-np.r_[0.5, 0.25] * spectral_width)
                ]
                .mean()
                .data.item()
                + detector_spectrum[
                    "t2" : tuple(np.r_[0.25, 0.5] * spectral_width)
                ]
                .mean()
                .data.item()
            ) / 2
            detector_spectra.append((detector_spectrum, raw_detector_spectrum))

        detector_spectrum, raw_detector_spectrum = detector_spectra[0]
        l2g_detector_spectrum, _ = detector_spectra[1]
        reference_maximum = detector_spectrum[
            "t2":inhomogeneous_bounds
        ].data.real.max()
        frequency_axis = detector_spectrum.getaxis("t2")
        inhomogeneous_center = inhomogeneous_bounds.mean()
        display_half_width = max(
            0.75 * inhomogeneous_linewidth,
            4 * homogeneous_linewidth,
        )

        fl.next(f"{experiment_type}: how both linewidths are determined")
        figure = plt.gcf()
        figure.clear()
        figure.set_size_inches(9, 9)
        detection_axis, comparison_axis = figure.subplots(
            2,
            1,
        )
        detection_axis.plot(
            frequency_axis,
            detector_spectrum.data.real / reference_maximum,
            color="k",
            label="conditioned, smoothed, baselined magnitude",
        )
        detection_axis.plot(
            frequency_axis,
            l2g_detector_spectrum.data.real
            / l2g_detector_spectrum.data.real.max(),
            color="tab:blue",
            alpha=0.65,
            label="equal-energy L2G coarse locator",
        )
        for threshold, linestyle in [(0.50, "--"), (0.10, "-."), (0.03, ":")]:
            detection_axis.axhline(
                threshold,
                color="tab:gray",
                linestyle=linestyle,
                label=f"{threshold:.0%} threshold",
            )
        for bound in inhomogeneous_bounds:
            detection_axis.axvline(bound, color="tab:red", linestyle="--")
        detection_axis.set_ylabel("normalized detector magnitude")
        detection_axis.set_title(
            "L2G locator: peak detail within +/-1200 Hz\n"
            "Complete final bounds are shown below"
        )
        detection_axis.set_xlabel("frequency / Hz")
        # Show the L2G peak-selection behavior over a physically relevant
        # region without allowing the full spectral width to flatten the peak.
        detection_axis.set_xlim(
            inhomogeneous_center - 1200,
            inhomogeneous_center + 1200,
        )
        detection_axis.grid(True)
        detection_axis.legend()

        raw_magnitude = raw_detector_spectrum.data.real.copy()
        raw_magnitude /= raw_magnitude.max()
        peak_frequency = raw_detector_spectrum.getaxis("t2")[
            raw_magnitude.argmax()
        ]
        single_lorentzian = 1 / (
            1
            + (
                2
                * (raw_detector_spectrum.getaxis("t2") - peak_frequency)
                / homogeneous_linewidth
            )
            ** 2
        )
        comparison_axis.plot(
            raw_detector_spectrum.getaxis("t2"),
            raw_magnitude,
            color="goldenrod",
            label="measured ensemble magnitude",
        )
        comparison_axis.plot(
            raw_detector_spectrum.getaxis("t2"),
            single_lorentzian,
            color="tab:blue",
            label="modeled absorptive Lorentzian (unit height)",
        )
        for bound in inhomogeneous_bounds:
            comparison_axis.axvline(
                bound,
                color="tab:red",
                linestyle="--",
                label=(
                    "inhomogeneous bounds"
                    if bound == inhomogeneous_bounds[0]
                    else None
                ),
            )
        for bound in peak_frequency + np.r_[-0.5, 0.5] * homogeneous_linewidth:
            comparison_axis.axvline(
                bound,
                color="tab:blue",
                linestyle=":",
                label=(
                    "homogeneous FWHM"
                    if bound == peak_frequency - 0.5 * homogeneous_linewidth
                    else None
                ),
            )
        comparison_axis.set_xlabel("frequency / Hz")
        comparison_axis.set_ylabel("normalized magnitude")
        comparison_axis.set_title(
            f"Envelope estimate: hom FWHM {homogeneous_linewidth:.1f} Hz\n"
            f"Ensemble extent at 3%: {inhomogeneous_linewidth:.1f} Hz"
        )
        comparison_axis.set_xlim(
            inhomogeneous_center - display_half_width,
            inhomogeneous_center + display_half_width,
        )
        comparison_axis.grid(True)
        comparison_axis.legend()
        figure.tight_layout()

        print(
            experiment_type,
            "echo center:",
            echo_center,
            "homogeneous linewidth:",
            homogeneous_linewidth,
            "provisional LSQ:",
            used_lsq_fallback,
            "inhomogeneous bounds:",
            inhomogeneous_bounds,
            "inhomogeneous linewidth:",
            inhomogeneous_linewidth,
        )
