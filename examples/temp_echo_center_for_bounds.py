r"""
Frequency-Bound Detection Requires the Correct Echo Center
==========================================================

For the same experimental FIR or ODNP-enhancement DCCT acquisition, compare
frequency-bound detection after exponential-correlation centering with a
modestly late time origin.  The late origin removes part of the high-SNR
beginning of the decay, so its Fourier transform gives biased final
inhomogeneous bounds.  Both rows contain the same acquired peak; only the
assumed location of zero time changes.

This is a temporary diagnostic example and is not intended for committing.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

from pyspecProcScripts import (
    det_inh_bounds,
    fid_side_from_echo,
    find_exponential_echo_center,
    fit_envelope,
    select_pathway,
)
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
        preliminary_bounds = determine_preliminary_bounds(
            raw_data, configuration
        )
        good_center = find_exponential_echo_center(
            sign_corrected(
                raw_data,
                configuration,
                preliminary_bounds["inh_bounds"],
            ),
            decay_rate=250,
        )
        bad_center = good_center + 0.5e-3
        homogeneous_linewidth = fit_envelope(
            select_pathway(
                fid_side_from_echo(
                    sign_corrected(
                        raw_data,
                        configuration,
                        preliminary_bounds["inh_bounds"],
                    ),
                    good_center,
                ),
                configuration["signal_pathway"],
            ),
            mult_two=True,
        )
        center_results = []
        # Fix sign correction and the L2G width across both trials to isolate
        # the effect of time origin.  This is a controlled sensitivity test,
        # not a measurement of the actual acquisition's centering error.
        centered_input = sign_corrected(
            raw_data, configuration, preliminary_bounds["inh_bounds"]
        )

        # Both rows reuse the same acquired data.  Only the time origin used
        # to construct the FID side changes.
        for label, echo_center in [
            ("good exponential-correlation center", good_center),
            ("perturbed center (+0.5 ms)", bad_center),
        ]:
            detector_input = centered_input.C
            detector_input.set_prop("echo_center", echo_center)
            det_inh_bounds(
                detector_input,
                0.10,
                peak_lowest_thresh=0.03,
                signal_pathway=configuration["signal_pathway"],
                apodization_linewidth=homogeneous_linewidth,
            )
            bounds = detector_input.get_prop("inh_bounds").copy()
            bound_description = f"{bounds[0]:.0f} to {bounds[1]:.0f} Hz"

            fid_spectrum = fid_side_from_echo(
                centered_input,
                echo_center,
            ).ft("t2")
            center_results.append(
                (label, echo_center, bounds, bound_description, fid_spectrum)
            )

        fl.next(f"{experiment_type}: t=0 changes the DCCT spectrum")
        figure = plt.gcf()
        figure.set_size_inches(10, 12)
        grid = figure.add_gridspec(2, 1, hspace=0.3)
        # Include both sets of final bounds; the old fixed view hid the ODNP
        # left edge and prevented visual comparison of the linewidths.
        common_plot_bounds = (
            min(result[2][0] for result in center_results) - 200,
            max(result[2][1] for result in center_results) + 200,
        )
        common_scale = max(
            abs(result[4]["t2":common_plot_bounds]).data.max()
            for result in center_results
        )
        for row, result in enumerate(center_results):
            label, echo_center, _, bound_description, fid_spectrum = result
            psd.DCCT(
                fid_spectrum["t2":common_plot_bounds],
                fig=figure,
                bbox=grid[row],
                custom_scaling=True,
                scaling_factor=common_scale,
                title=(
                    f"{label}: t=0 at {echo_center * 1e3:.3f} ms; "
                    f"bounds {bound_description}"
                ),
            )

        fl.next(f"{experiment_type}: good vs bad frequency bounds")
        figure = plt.gcf()
        figure.clear()
        figure.set_size_inches(9, 7)
        axes = figure.subplots(2, 1, sharex=True, sharey=True)
        figure.suptitle(experiment_type)
        for axis, result in zip(axes, center_results):
            label, echo_center, bounds, bound_description, fid_spectrum = (
                result
            )
            signal_spectrum = select_pathway(
                fid_spectrum,
                configuration["signal_pathway"],
            ).mean(configuration["indirect"])
            axis.plot(
                signal_spectrum.getaxis("t2"),
                abs(signal_spectrum.data),
                color="k",
            )
            for bound in bounds:
                axis.axvline(bound, color="r", linestyle="--")
            axis.set_xlim(common_plot_bounds)
            axis.set_title(
                f"{label}: {echo_center * 1e3:.3f} ms\nfinal 3% bounds "
                f"{bound_description}"
            )
            axis.set_ylabel("signal magnitude")
            axis.grid(True)
            print(
                experiment_type,
                label,
                "3% bounds / Hz:",
                bounds,
                "width / Hz:",
                np.diff(bounds).item(),
            )
        axes[-1].set_xlabel("frequency / Hz")
        figure.tight_layout()

        print(
            experiment_type,
            "exponential-correlation center:",
            good_center,
            "bad center:",
            bad_center,
        )
