r"""
Spectral Wrapping When Alignment Follows a Narrow Frequency Slice
===============================================================

Use the FIR and ODNP data and known imposed transient shifts from the other
temporary examples.  Compare shifting the full spectrum before cropping with
cropping first and shifting on the resulting periodic Fourier grid.  Show
the untouched experimental DCCT data and the shifted input as references.

The diagnostic also shifts an identically masked spectrum on the original
wide grid.  It retains the same input samples as the narrow calculation but
lets outgoing signal leave the window.  Their difference isolates folding
at the narrow window boundary from missing signal due to initial clipping.
No estimated correlation corrections are used in this comparison.  The known
corrections are rounded to the nearest frequency bin (at most 1 Hz here), so
the boundary transfer can be verified without fractional-bin interpolation.

This temporary example is not intended for committing.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

from temp_bounds_dcct_data import (
    EXPERIMENT_TYPES,
    determine_demo_bounds,
    load_dcct_dataset,
    show_raw_data,
)

# {{{ changeable parameters
manual_half_widths = {
    "FIR": 400.0,
    "FIR no power": 400.0,
    "ODNP enhancement": 500.0,
}  # Hz
# }}}

with psd.figlist_var() as fl:
    for experiment_type in EXPERIMENT_TYPES:
        raw_data, shifted_data, configuration = load_dcct_dataset(
            experiment_type
        )
        show_raw_data(fl, raw_data, configuration)
        fl.next(f"{experiment_type}: imposed drift before correction")
        psd.DCCT(shifted_data, fig=plt.gcf(), title="imposed transient shifts")
        bounds = (
            determine_demo_bounds(raw_data, configuration)["inh_bounds"].mean()
            + np.r_[-1, 1] * manual_half_widths[experiment_type]
        )
        narrow_input = shifted_data["t2":bounds].C
        # Match the mask to actual retained bins, not nominal slice endpoints.
        retained_bounds = narrow_input.getaxis("t2")[[0, -1]]
        masked_input = shifted_data * shifted_data.fromaxis("t2").run(
            lambda frequency: (frequency >= retained_bounds[0])
            & (frequency <= retained_bounds[1])
        )
        frequency_bin = abs(np.diff(narrow_input.getaxis("t2")[:2]).item())
        known_corrections = configuration["frequency_shifts"].C
        known_corrections.run(
            lambda shift: np.rint(shift / frequency_bin) * frequency_bin
        )
        corrected_spectra = []
        for input_spectrum in [shifted_data, narrow_input, masked_input]:
            corrected = input_spectrum.C
            corrected.set_error(None)
            corrected.ift("t2").ift(configuration["phase_dimensions"])
            corrected *= np.exp(
                -1j * 2 * np.pi * known_corrections * corrected.fromaxis("t2")
            )
            corrected.ft("t2").ft(configuration["phase_dimensions"])
            corrected_spectra.append(corrected)
        full_reference = corrected_spectra[0]["t2":retained_bounds].C
        wrapped_spectrum = corrected_spectra[1]
        masked_reference = corrected_spectra[2]["t2":retained_bounds].C

        # Inspect one transient with the largest outgoing energy.  Choose it
        # by a stated numerical criterion so the plotted edge transfer is
        # reproducible and is not mistaken for an ensemble-average trace.
        transient_dimensions = [
            dim for dim in shifted_data.dimlabels if dim != "t2"
        ] + ["t2"]
        wide_transients = (
            corrected_spectra[2]
            .C.ift(configuration["phase_dimensions"])
            .reorder(transient_dimensions)
        )
        narrow_transients = wrapped_spectrum.C.ift(
            configuration["phase_dimensions"]
        ).reorder(transient_dimensions)
        wide_frequency = wide_transients.getaxis("t2")
        wide_traces = wide_transients.data.reshape(-1, wide_frequency.size)
        narrow_traces = narrow_transients.data.reshape(
            -1, narrow_transients.shape["t2"]
        )
        outside_window = (wide_frequency < retained_bounds[0]) | (
            wide_frequency > retained_bounds[1]
        )
        trace_index = np.argmax(
            np.sum(abs(wide_traces[:, outside_window]) ** 2, axis=1)
        )
        outgoing_index = np.argmax(
            abs(wide_traces[trace_index]) * outside_window
        )
        outgoing_frequency = wide_frequency[outgoing_index]
        window_period = narrow_transients.shape["t2"] * frequency_bin
        folded_frequency = (
            outgoing_frequency - retained_bounds[0]
        ) % window_period + retained_bounds[0]
        # Explicitly fold every wide-grid bin into the narrow period.  This
        # independent bin-summing check verifies actual wrapping, rather than
        # inferring aliasing merely from a large correlation error.
        folded_trace = np.zeros(narrow_transients.shape["t2"], dtype=complex)
        # A nonzero recorded time origin gives each spectral replica a phase
        # factor.  Preserve that metadata in the independent complex sum.
        fold_number = np.floor(
            (wide_frequency - retained_bounds[0]) / window_period
        )
        np.add.at(
            folded_trace,
            np.rint(
                (wide_frequency - retained_bounds[0]) / frequency_bin
            ).astype(int)
            % folded_trace.size,
            wide_traces[trace_index]
            * np.exp(
                1j
                * 2
                * np.pi
                * fold_number
                * window_period
                * narrow_input.C.ift("t2").getaxis("t2")[0]
            ),
        )
        np.testing.assert_allclose(
            narrow_traces[trace_index],
            folded_trace,
            atol=1e-8 * abs(folded_trace).max(),
            rtol=1e-8,
        )
        fl.next(f"{experiment_type}: signal crosses the frequency boundary")
        figure = plt.gcf()
        figure.set_size_inches(9, 7)
        outgoing_axis, folded_axis = figure.subplots(
            2, 1, sharex=True, sharey=True
        )
        outgoing_axis.plot(
            wide_frequency,
            abs(wide_traces[trace_index]),
            color="k",
            label="same masked input corrected on wide grid",
        )
        folded_axis.plot(
            narrow_transients.getaxis("t2"),
            abs(narrow_traces[trace_index]),
            color="tab:red",
            label="correction on narrow grid",
        )
        folded_axis.plot(
            narrow_transients.getaxis("t2"),
            abs(folded_trace),
            "--",
            color="k",
            label="independent sum of folded wide-grid bins",
        )
        outgoing_axis.axvline(outgoing_frequency, color="tab:blue", ls=":")
        folded_axis.axvline(folded_frequency, color="tab:blue", ls=":")
        outgoing_axis.set_title(
            f"{experiment_type}: transient {trace_index}, "
            "largest outgoing energy"
        )
        folded_axis.set_title(
            f"Outgoing signal at {outgoing_frequency:.0f} Hz folds to "
            f"{folded_frequency:.0f} Hz (period {window_period:.0f} Hz)"
        )
        for axis in [outgoing_axis, folded_axis]:
            for bound in retained_bounds:
                axis.axvline(bound, color="tab:gray", ls="--")
            axis.set_xlim(
                retained_bounds[0] - configuration["max_alignment_shift"],
                retained_bounds[1] + configuration["max_alignment_shift"],
            )
            axis.set_ylabel("transient magnitude")
            axis.legend(fontsize=9)
            axis.grid(True)
        folded_axis.set_xlabel("frequency / Hz")
        figure.tight_layout()

        fl.next(f"{experiment_type}: full correction vs narrow-grid wrapping")
        figure = plt.gcf()
        figure.set_size_inches(10, 12)
        grid = figure.add_gridspec(2, 1, hspace=0.3)
        common_scale = max(
            abs(spectrum).data.max()
            for spectrum in [full_reference, wrapped_spectrum]
        )
        for row, (label, spectrum) in enumerate(
            [
                ("correct on full spectrum, then crop", full_reference),
                (
                    "crop first, then correct: periodic wrapping",
                    wrapped_spectrum,
                ),
            ]
        ):
            psd.DCCT(
                spectrum,
                fig=figure,
                bbox=grid[row],
                custom_scaling=True,
                scaling_factor=common_scale,
                title=f"{label}\n{retained_bounds[0]:.0f} to "
                f"{retained_bounds[1]:.0f} Hz; known shifts",
            )

        # Compare magnitudes only after forming the complex difference.  This
        # preserves contributions whose phases cancel in magnitude subtraction.
        folded_component = wrapped_spectrum - masked_reference
        fl.next(f"{experiment_type}: isolate the folded contribution")
        figure = plt.gcf()
        figure.set_size_inches(9, 7)
        comparison_axis, difference_axis = figure.subplots(2, 1, sharex=True)
        for label, spectrum in [
            ("full-spectrum correction", full_reference),
            ("same clipped input, shifted on wide grid", masked_reference),
            ("same clipped input, shifted on narrow grid", wrapped_spectrum),
        ]:
            envelope = abs(spectrum).mean_all_but(["t2"])
            comparison_axis.plot(
                envelope.getaxis("t2"), envelope.data.real, label=label
            )
        folded_envelope = abs(folded_component).mean_all_but(["t2"])
        difference_axis.plot(
            folded_envelope.getaxis("t2"),
            folded_envelope.data.real,
            color="tab:red",
            label="magnitude of narrow-grid minus wide-grid correction",
        )
        for axis in [comparison_axis, difference_axis]:
            axis.set_xlim(retained_bounds)
            axis.set_ylabel("mean magnitude over pathways")
            axis.grid(True)
            axis.legend(fontsize=9)
        comparison_axis.set_title(
            f"{experiment_type}: identical known corrections in all paths"
        )
        difference_axis.set_xlabel("frequency / Hz")
        figure.tight_layout()
        print(
            experiment_type,
            "window / Hz:",
            retained_bounds,
            "folded-component RMS / full-reference RMS:",
            np.linalg.norm(folded_component.data)
            / np.linalg.norm(full_reference.data),
        )
