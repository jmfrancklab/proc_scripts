r"""
Frequency Bounds Control DCCT Correlation Alignment
===================================================

For the same experimental FIR or ODNP-enhancement acquisition, compare
manual +/-400, +/-450, and +/-500 Hz windows with linewidth-aware bounds.
All paths receive identical data, imposed shifts, zeroth-order phasing, and
sign correction.  A shared phase is estimated from the signal pathway at the
unshifted echo center before any frequency slicing.  The complete echo is
retained for correlation, and signed FIR/enhancement amplitudes are preserved
in the displayed data.
The manual window is centered at the midpoint of the detected inhomogeneous
bounds.  This is the actual processing slice, not just the displayed x range.

All candidate widths are shown because the alignment error need not change
monotonically with window width.  The recovered corrections are applied to
the full spectrum here; temp_frequency_window_wrapping.py separately shows
spectral wrapping when corrections are applied inside a cropped window.

This is a temporary diagnostic example and is not intended for committing.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

from pyspecProcScripts import (
    correl_align,
    fid_side_from_echo,
    select_pathway,
    zeroth_order_ph,
)
from temp_bounds_dcct_data import (
    EXPERIMENT_TYPES,
    determine_demo_bounds,
    frequency_mask,
    keep_all_coherence_pathways,
    load_dcct_dataset,
    show_raw_data,
    sign_corrected,
)

# {{{ changeable parameters
manual_half_widths = (400.0, 450.0, 500.0)  # Hz
# }}}

plt.rcParams["image.aspect"] = "auto"
with psd.figlist_var() as fl:
    for experiment_type in EXPERIMENT_TYPES:
        raw_data, data, configuration = load_dcct_dataset(experiment_type)
        show_raw_data(fl, raw_data, configuration)
        detected_bounds = determine_demo_bounds(raw_data, configuration)
        # Estimate one phase from the reference acquisition at the interpolated
        # echo center.  Only the reference copy is sliced to obtain that point;
        # the compared spectra retain their full echoes.  A common rotation
        # preserves the phase-cycle relationships and cancels in correlation.
        phase_correction = zeroth_order_ph(
            select_pathway(
                fid_side_from_echo(raw_data, detected_bounds["echo_center"]),
                configuration["signal_pathway"],
            )["t2", 0]
        )
        data /= phase_correction
        print(
            experiment_type,
            "removed zeroth-order phase / degrees:",
            np.angle(phase_correction, deg=True),
        )
        fl.next(f"{experiment_type}: imposed drift before alignment")
        psd.DCCT(
            data,
            fig=plt.gcf(),
            title="zeroth-order phased signal with imposed transient shifts",
        )
        # Determine signs once before slicing, so the comparison changes only
        # the processing window and does not also change sign estimation.
        alignment_data = sign_corrected(
            data, configuration, detected_bounds["inh_bounds"]
        )
        alignment_results = []

        # Both paths start from copies of the exact same DCCT dataset.  The
        # only difference is whether the slice includes the homogeneous tails
        # and the maximum expected alignment shift.
        for label, bounds in [
            (
                "linewidth-aware bounds",
                detected_bounds["processing_bounds"],
            ),
        ] + [
            (
                f"manual ±{half_width:.0f} Hz bounds",
                detected_bounds["inh_bounds"].mean()
                + np.r_[-1, 1] * half_width,
            )
            for half_width in manual_half_widths
        ]:
            alignment_input = alignment_data["t2":bounds].C
            estimated_shifts = correl_align(
                alignment_input,
                frq_mask_fn=frequency_mask(bounds),
                coherence_unmask_fn=keep_all_coherence_pathways,
                repeat_dims=[configuration["indirect"]],
                max_shift=configuration["max_alignment_shift"],
                fig_title=f"{experiment_type}: {label}",
            )
            # The correlation determines relative shifts, so remove the
            # arbitrary common frequency offset before comparing with truth.
            estimated_shifts -= estimated_shifts.data.mean().item()
            aligned = data.C.ift("t2").ift(configuration["phase_dimensions"])
            aligned *= np.exp(
                -1j * 2 * np.pi * estimated_shifts * aligned.fromaxis("t2")
            )
            aligned.ft("t2").ft(configuration["phase_dimensions"])
            alignment_results.append(
                (label, bounds, estimated_shifts, aligned)
            )

        common_plot_bounds = (-1800.0, 2000.0)
        common_scale = max(
            abs(aligned["t2":common_plot_bounds]).data.max()
            for _, _, _, aligned in alignment_results
        )
        # Keep readable top/bottom maps for every candidate, with the same
        # broad reference and intensity scale across all comparisons.
        for candidate in alignment_results[1:]:
            fl.next(f"{experiment_type}: broad vs {candidate[0]} DCCT")
            figure = plt.gcf()
            figure.set_size_inches(10, 12)
            grid = figure.add_gridspec(2, 1, hspace=0.3)
            for row, (label, bounds, _, aligned) in enumerate(
                [alignment_results[0], candidate]
            ):
                psd.DCCT(
                    aligned["t2":common_plot_bounds],
                    fig=figure,
                    bbox=grid[row],
                    custom_scaling=True,
                    scaling_factor=common_scale,
                    title=f"{label}: {bounds[0]:.0f} to {bounds[1]:.0f} Hz",
                )

        true_shifts = configuration["frequency_shifts"].data.copy()
        true_shifts -= true_shifts.mean()
        fl.next(f"{experiment_type}: recovered frequency shifts")
        plt.plot(
            true_shifts.ravel(),
            true_shifts.ravel(),
            "k--",
            label="ideal",
        )
        for label, _, estimated_shifts, _ in alignment_results:
            shift_error = estimated_shifts.data - true_shifts
            print(
                experiment_type,
                label,
                "RMS shift error / Hz:",
                np.sqrt(np.mean(shift_error**2)),
            )
            plt.plot(
                true_shifts.ravel(),
                estimated_shifts.data.ravel(),
                "o",
                label=(
                    f"{label}; RMS error "
                    f"{np.sqrt(np.mean(shift_error**2)):.1f} Hz"
                ),
            )
        plt.xlabel("imposed relative shift / Hz")
        plt.ylabel("correlation-alignment result / Hz")
        plt.legend()
        plt.grid(True)
