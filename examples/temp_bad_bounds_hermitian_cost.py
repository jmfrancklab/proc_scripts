r"""
Frequency Bounds Control the Hermitian Cost Function
====================================================

For the same experimental FIR or ODNP-enhancement DCCT acquisition, compare the
Hermitian cost obtained from the linewidth-aware processing bounds with the
cost obtained from a practical manually chosen ±1800 Hz window.  A rectangular
frequency window corresponds to convolution with a sinc-like kernel in time.
The narrow slice also lengthens the dwell and the three-dwell exclusion at
the start of the Hermitian search.  Both effects can displace its minimum;
this comparison includes the complete processing effect of choosing bounds.
Unlike the earlier half-height example, both windows retain an interior
minimum.  These measured data need not exhibit a catastrophic failure: the
comparison shows the actual change in the minimum and its surrounding shape.
Both estimates use the same sign-corrected average and acquisition-time origin.

This is a temporary diagnostic example and is not intended for committing.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

from pyspecProcScripts import (
    find_exponential_echo_center,
    hermitian_function_test,
    select_pathway,
)
from temp_bounds_dcct_data import (
    EXPERIMENT_TYPES,
    determine_demo_bounds,
    load_dcct_dataset,
    show_raw_data,
    sign_corrected,
)

plt.rcParams["image.aspect"] = "auto"
with psd.figlist_var() as fl:
    for experiment_type in EXPERIMENT_TYPES:
        raw_data, _, configuration = load_dcct_dataset(experiment_type)
        show_raw_data(fl, raw_data, configuration)
        detected_bounds = determine_demo_bounds(raw_data, configuration)
        hermitian_results = []
        averaged_signal = select_pathway(
            sign_corrected(
                raw_data, configuration, detected_bounds["inh_bounds"]
            ),
            configuration["signal_pathway"],
        ).mean(configuration["indirect"])
        # FIR/enhancement sign changes alter a complex average.  Compare the
        # two estimators on this same average, not on differently signed data.
        reference_center = find_exponential_echo_center(
            averaged_signal, decay_rate=250
        )
        if averaged_signal.get_prop("dig_filter") is not None:
            averaged_signal *= averaged_signal.get_prop("dig_filter")
        acquisition_start = averaged_signal.C.ift("t2").getaxis("t2")[0]

        fl.next(f"{experiment_type}: DCCT inputs for Hermitian phasing")
        figure = plt.gcf()
        figure.set_size_inches(10, 12)
        grid = figure.add_gridspec(2, 1, hspace=0.3)
        for row, (label, bounds) in enumerate(
            [
                (
                    "good linewidth-aware bounds",
                    detected_bounds["processing_bounds"],
                ),
                (
                    "manual ±1800 Hz bounds",
                    tuple(
                        np.mean(detected_bounds["inh_bounds"])
                        + np.array([-1800, 1800])
                    ),
                ),
            ]
        ):
            psd.DCCT(
                raw_data["t2":bounds],
                fig=figure,
                bbox=grid[row],
                title=f"{label}: {bounds[0]:.0f} to {bounds[1]:.0f} Hz",
            )

            hermitian_input = averaged_signal["t2":bounds].C.ift("t2")
            # Cropping changes the time grid, including its first coordinate.
            # Convert the returned center back to the full acquisition origin.
            origin_offset = (
                hermitian_input.getaxis("t2")[0] - acquisition_start
            )
            hermitian_input["t2"] -= hermitian_input.getaxis("t2")[0]
            aliasing_slop = 3
            # Keep the routine's mixed energy/envelope diagnostics separate
            # from the presentation figures; extract its actual cost, without
            # reimplementing or smoothing the calculation.
            cost_figures = psd.figlist_var()
            # An independent figlist otherwise starts numbering at figure 1,
            # which would reuse the raw-data figure in the outer figlist.
            cost_figures.next("power terms", fig=plt.figure())
            detected_center = (
                hermitian_function_test(
                    hermitian_input,
                    aliasing_slop=aliasing_slop,
                    fl=cost_figures,
                )
                + origin_offset
            )
            # Read this call's diagnostic directly instead of scanning every
            # open figure (which can accidentally select an earlier curve).
            cost_figures.next("power terms")
            cost_line = next(
                line
                for line in plt.gca().lines
                if line.get_label() == "cost function"
            )
            cost_time_to_ms = {
                "s": 1e3,
                "ms": 1,
                r"\mu s": 1e-3,
            }[cost_figures.units[cost_figures.current]]
            cost_values = np.asarray(cost_line.get_ydata()).copy()
            assert np.all(np.isfinite(cost_values))
            assert 0 < np.argmin(cost_values) < cost_values.size - 1
            # The diagnostic axis excludes the rounded three-dwell offset.
            # Recover that exact offset from the returned minimum rather than
            # assuming the oversampled rounding leaves it unchanged.
            cost_times = np.asarray(cost_line.get_xdata()).copy()
            cost_times *= cost_time_to_ms
            cost_times += (
                detected_center * 1e3 - cost_times[np.argmin(cost_values)]
            )
            assert np.isclose(
                cost_times[np.argmin(cost_values)], detected_center * 1e3
            )
            hermitian_results.append(
                (
                    label,
                    bounds,
                    detected_center,
                    cost_times,
                    cost_values,
                    cost_times[0],
                )
            )
            plt.close(plt.gcf())

        fl.basename = None
        fl.next(f"{experiment_type}: Hermitian cost comparison")
        figure = plt.gcf()
        figure.clear()
        figure.set_size_inches(9, 8)
        axes = figure.subplots(2, 1, sharex=True, sharey=True)
        figure.suptitle(f"{experiment_type}: Hermitian cost near the echo")
        for axis, result in zip(axes, hermitian_results):
            label, bounds, detected_center, cost_time, cost, search_start = (
                result
            )
            cost /= cost.max()
            axis.plot(
                cost_time, cost, color="tab:purple", label="measured cost"
            )
            axis.axvline(
                reference_center * 1e3,
                color="k",
                linestyle="--",
                label="exponential-correlation estimate",
            )
            axis.axvline(
                detected_center * 1e3,
                color="r",
                linestyle=":",
                label=f"minimum: {detected_center * 1e3:.3f} ms",
            )
            axis.axvline(
                search_start,
                color="tab:gray",
                linestyle="-.",
                label="start of searched cost interval",
            )
            axis.plot(detected_center * 1e3, cost.min(), "o", color="tab:red")
            axis.set_title(f"{label}: {bounds[0]:.0f} to {bounds[1]:.0f} Hz")
            axis.set_ylabel("cost / maximum over search")
            axis.legend(loc="upper right", fontsize=9)
            axis.set_xlim(0.4, 2.5)
            axis.set_ylim(bottom=0)
            axis.grid(True)
            print(
                experiment_type,
                label,
                "Hermitian center / ms:",
                detected_center * 1e3,
                "searched interval starts / ms:",
                search_start,
            )
        axes[-1].set_xlabel(
            "candidate echo center since acquisition start / ms"
        )
        figure.tight_layout()
