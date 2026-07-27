from ..calc_error import calc_masked_variance
from ..correlation_alignment import correl_align
from ..phasing import (
    determine_sign,
    fid_from_echo,
    find_peakrange,
    hermitian_function_test,
    zeroth_order_ph,
)
from ..simple_functions import select_pathway
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, r_
import pyspecdata as psd


def table_of_integrals(
    s,
    signal_range=None,
    signal_pathway=None,
    fl=None,
    echo_like=True,
    title="",
    direct="t2",
    expansion=2,
    peak_lower_thresh=0.1,
    inc_plot_color=True,
    repeat_dims=None,
    non_repeat_dims=None,
    alignment_sigma=150.0,
    max_shift_Hz=880.0,
    center_aligned_peak=True,
    show_alignment_diagnostics=False,
    propagate_error=True,
    excluded_pathways=None,
    fallback_signal=None,
    fid_from_echo_slice_multiplier=5,
):
    """Generate a table of frequency-domain integrals after DCCT alignment.

    The echo-like path follows the FIR DCCT pipeline used for the ODNP T1
    examples: remove receiver offset in the time domain, Hermitian-center the
    echo, use an exponential echo filter only for correlation alignment, slice
    the FID from the aligned echo, and integrate the selected coherence
    pathway.  Setting ``echo_like=False`` preserves the older direct-spectrum
    behavior while still using correlation alignment.
    """

    def mean_if_present(x, dimnames):
        for dimname in dimnames:
            if dimname in x.dimlabels:
                x = x.mean(dimname)
        return x

    def to_time_domain(x):
        out = x.C
        if out.get_ft_prop(direct):
            out.ift(direct)
        return out

    def to_frequency_domain(x):
        out = x.C
        if not out.get_ft_prop(direct):
            out.ft(direct)
        return out

    if fl is None:
        raise ValueError("table_of_integrals requires a figlist via fl")
    signal_pathway = (
        s.get_prop("coherence_pathway")
        if signal_pathway is None
        else signal_pathway
    )
    if signal_pathway is None:
        raise ValueError(
            "table_of_integrals needs signal_pathway or a coherence_pathway"
            " property"
        )
    s = s.C
    s.set_prop("coherence_pathway", signal_pathway)
    if non_repeat_dims is None:
        non_repeat_dims = []
    elif isinstance(non_repeat_dims, str):
        non_repeat_dims = [non_repeat_dims]
    else:
        non_repeat_dims = list(non_repeat_dims)
    if repeat_dims is None:
        phcycdims = [j for j in s.dimlabels if j.startswith("ph")]
        repeat_dims = [
            j
            for j in s.dimlabels
            if j not in set(phcycdims + [direct, "nScans"])
            and j not in non_repeat_dims
        ]
    elif isinstance(repeat_dims, str):
        repeat_dims = [repeat_dims]
    else:
        repeat_dims = list(repeat_dims)
    if len(repeat_dims) == 0:
        raise ValueError(
            "table_of_integrals needs at least one repeat dimension for"
            " correlation alignment"
        )

    # {{{ Determine the initial signal range
    # The fallback signal is used only to get through Hermitian centering,
    # alignment, and FID slicing when the current low-SNR node does not pass
    # find_peakrange.  For echo-like data, the final integration range is
    # recalculated from the filtered/aligned current node below.
    if signal_range is None:
        signal_for_range = (
            fallback_signal
            if fallback_signal is not None
            else select_pathway(s.C, signal_pathway)
        )
        frq_center, frq_half = find_peakrange(
            signal_for_range,
            fl=None,
            direct=direct,
            peak_lower_thresh=peak_lower_thresh,
        )
        signal_range = tuple(sorted(frq_center + r_[-1, 1] * abs(frq_half)))
    else:
        if signal_range == "peakrange":
            signal_range = s.get_prop("peakrange")
        frq_center = np.mean(signal_range)
        frq_half = abs(0.5 * np.diff(signal_range).item())
    signal_range_expanded = tuple(
        sorted(frq_center + expansion * r_[-1, 1] * frq_half)
    )
    # }}}

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.suptitle(title)
    fl.next("Raw Data with correlation alignment", fig=fig)
    fl.skip_units_check()

    # {{{ Receiver-offset correction and Hermitian centering
    working = to_time_domain(s)
    working.set_units(direct, "s")
    if echo_like:
        working.ift(list(signal_pathway))
        t_max = working.getaxis(direct)[-1]
        working -= working[direct : (t_max * 0.75, None)].mean(direct)
        working.ft(list(signal_pathway))
    working.setaxis(direct, lambda x: x - x[0])
    hermitian_input = select_pathway(working.C, signal_pathway)
    hermitian_input = mean_if_present(
        hermitian_input, repeat_dims + ["nScans", "repeats"]
    )
    hermitian_shift = hermitian_function_test(
        hermitian_input,
        direct=direct,
        fl=fl if show_alignment_diagnostics else None,
    )
    working.setaxis(direct, lambda x: x - hermitian_shift)
    working.register_axis({direct: 0})
    working /= zeroth_order_ph(
        select_pathway(working[direct:0], signal_pathway)
    )
    # }}}

    if echo_like:
        fid_unapodized = working.C[direct : (0, None)]
        fid_unapodized *= 2
        fid_unapodized[direct:0] *= 0.5
        fid_unapodized.ft(direct)
        alignment_data = to_time_domain(fid_unapodized)
        acq_params = alignment_data.get_prop("acq_params")
        actual_tau = acq_params["tau_us"] * 1e-6
        alignment_data *= np.exp(
            -abs(alignment_data.fromaxis(direct) - actual_tau) / 10e-3
        )
        alignment_data.ft(direct)
        if fallback_signal is not None:
            # {{{ Recalculate current-node range after alignment filtering
            # The fallback range is allowed to
            # reach the alignment/FID-slicing stage, but the final limits are
            # derived from the current node after the exponential echo filter.
            alignment_signal_for_integral = select_pathway(
                alignment_data.C, signal_pathway
            )
            alignment_signal_for_integral = mean_if_present(
                alignment_signal_for_integral, ("nScans", "repeats")
            )
            if len(repeat_dims) == 1:
                alignment_signal_for_integral = alignment_signal_for_integral[
                    repeat_dims[0], -1
                ]
            argmax_frq = (
                alignment_signal_for_integral.C.run(abs).argmax(direct).item()
            )
            peak_search_slice = tuple(
                sorted(argmax_frq + r_[-1, 1] * abs(frq_half) / 2)
            )
            frq_center, frq_half = find_peakrange(
                alignment_signal_for_integral[direct:peak_search_slice],
                direct=direct,
                peak_lower_thresh=peak_lower_thresh,
                fl=None,
            )
            frq_half = abs(frq_half)
            signal_range = tuple(sorted(frq_center + r_[-1, 1] * frq_half))
            signal_range_expanded = tuple(
                sorted(frq_center + expansion * r_[-1, 1] * frq_half)
            )
            print(
                "fallback peak range was only used to reach FID slicing; "
                "integration limits were recalculated from the filtered "
                "current node"
            )
            # }}}
    else:
        alignment_data = to_frequency_domain(working)

    # {{{ Correlation alignment
    def frq_mask(x):
        signal_pathway_for_mask = x.get_prop("coherence_pathway")
        signal = select_pathway(x, signal_pathway_for_mask)
        signal = mean_if_present(signal, ("nScans", "repeats"))
        nu_center = signal.argmax(direct)
        return x * np.exp(
            -((x.fromaxis(direct) - nu_center) ** 2) / (4 * alignment_sigma**2)
        )

    def coherence_unmask(coh_array):
        if len(signal_pathway) == 0:
            coh_array.data[:] = 1
            return coh_array
        for ph_name, ph_val in coh_array.get_prop("coherence_pathway").items():
            coh_array[ph_name, ph_val] = 1
        return coh_array

    mysign_for_alignment = (
        select_pathway(alignment_data[direct:signal_range].C, signal_pathway)
        .real.sum(direct)
        .run(np.sign)
    )
    mysign_for_alignment[lambda x: x == 0] = 1
    alignment_input = alignment_data.C
    alignment_input.reorder(direct, first=False)
    opt_shift = correl_align(
        alignment_input.C * mysign_for_alignment,
        frq_mask_fn=frq_mask,
        coherence_unmask_fn=coherence_unmask,
        repeat_dims=repeat_dims,
        non_repeat_dims=non_repeat_dims,
        max_shift=max_shift_Hz,
        direct=direct,
        fl=fl if show_alignment_diagnostics else None,
    )
    aligned = working.C
    aligned.ft(direct)
    aligned.ift(list(signal_pathway))
    aligned.ift(direct)
    aligned *= np.exp(-1j * 2 * pi * opt_shift * aligned.fromaxis(direct))
    aligned.ft(list(signal_pathway))
    aligned.ft(direct)
    # }}}

    if center_aligned_peak:
        aligned_peak = abs(
            select_pathway(
                aligned[direct:signal_range_expanded].C, signal_pathway
            )
        )
        aligned_peak.mean_all_but([direct])
        frq_center = aligned_peak.argmax(direct).item()
        signal_range = tuple(sorted(frq_center + r_[-1, 1] * frq_half))
        signal_range_expanded = tuple(
            sorted(frq_center + expansion * r_[-1, 1] * frq_half)
        )

    selected = select_pathway(
        aligned[direct:signal_range_expanded], signal_pathway
    )
    fl.image(selected, ax=ax1, interpolation="auto")
    ax1.set_title("extract signal pathway\nand correlation align")

    mysign = determine_sign(
        selected.C,
        signal_range,
        direct=direct,
    )
    selected *= mysign
    fl.image(selected, ax=ax2, interpolation="auto")
    ax2.set_title("check phase variation\nalong indirect")

    if echo_like:
        aligned_fid = fid_from_echo(
            aligned.C.set_error(None),
            signal_pathway,
            frq_center=frq_center,
            frq_half=frq_half,
            fl=fl if show_alignment_diagnostics else None,
            direct=direct,
            slice_multiplier=fid_from_echo_slice_multiplier,
        )
        aligned_fid = to_frequency_domain(aligned_fid)
        selected = select_pathway(
            aligned_fid[direct:signal_range_expanded], signal_pathway
        )
    else:
        selected *= mysign

    fl.image(selected, ax=ax3, interpolation="auto")
    ax3.set_title(
        "FID sliced" + (", phased," if echo_like else "") + " and aligned"
    )

    if echo_like:
        d = aligned_fid.C
        d.ift(direct)
        d /= zeroth_order_ph(select_pathway(d[direct:0], signal_pathway))
        d.ft(direct)
        for ph_name in signal_pathway:
            d.set_ft_prop(ph_name, "unitary", True)
        this_IR = mean_if_present(d.C, ("nScans",))
        this_IR.set_prop("coherence_pathway", signal_pathway)
        excluded_pathways_for_error = []
        if excluded_pathways is not None:
            excluded_pathways_for_error.extend(excluded_pathways)
        if set(["ph1", "ph2"]).issubset(signal_pathway):
            excluded_pathways_for_error.extend(
                [
                    {"ph1": 0, "ph2": 0},
                    {
                        "ph1": signal_pathway["ph1"],
                        "ph2": signal_pathway["ph2"],
                    },
                ]
            )
        else:
            excluded_pathways_for_error.append(signal_pathway)
        signal_for_integral = select_pathway(this_IR.C, signal_pathway)
        signal_for_integral = mean_if_present(
            signal_for_integral, repeat_dims + ["nScans", "repeats"]
        )
        frq_center, frq_half = find_peakrange(
            signal_for_integral,
            direct=direct,
            peak_lower_thresh=peak_lower_thresh,
            fl=None,
        )
        frq_half = abs(frq_half)
        peak_frq_slice = list(sorted(frq_center + r_[-1, 1] * frq_half))
        frq_slice = list(
            sorted(
                frq_center
                + r_[-1, 1] * frq_half * fid_from_echo_slice_multiplier
            )
        )
        df = this_IR.get_ft_prop(direct, "df")
        for j in range(2):
            idx = np.searchsorted(this_IR[direct], frq_slice[j] + 0.5 * df)
            idx = min(idx, len(this_IR[direct]) - 1)
            frq_slice[j] = this_IR[direct][idx]
        integral_error = None
        if propagate_error:
            spectral_datapoint_variance = calc_masked_variance(
                this_IR,
                excluded_frqs=[peak_frq_slice],
                indirect=repeat_dims,
                excluded_pathways=excluded_pathways_for_error,
                direct=direct,
            )
            points_in_slice = this_IR[direct:frq_slice].shape[direct]
            integral_error = np.sqrt(
                spectral_datapoint_variance.data * df**2 * points_in_slice
            )
        selected = select_pathway(
            this_IR[direct:frq_slice], signal_pathway
        ).integrate(direct)
        selected.set_error(integral_error)
        fl.next("IR Integration Limits")
        for j in range(len(d.getaxis(repeat_dims[0]))):
            fl.plot(
                select_pathway(
                    mean_if_present(d[repeat_dims[0], j].C, ("nScans",)),
                    signal_pathway,
                ),
                label=f"{j}",
                human_units=False,
            )
        plt.axvline(frq_slice[0])
        plt.axvline(frq_slice[-1])
    else:
        final_slice = selected[direct:signal_range]
        if final_slice.shape[direct] < 2:
            final_peak = abs(selected.C)
            final_peak.mean_all_but([direct])
            frq_center = final_peak.argmax(direct).item()
            signal_range = tuple(sorted(frq_center + r_[-1, 1] * frq_half))
            final_slice = selected[direct:signal_range]
            if final_slice.shape[direct] < 2:
                raise ValueError(
                    "table_of_integrals could not find enough points in the "
                    "final integration window"
                )
        selected = final_slice.real.integrate(direct).set_error(None)
    if inc_plot_color:
        selected.set_plot_color_next()
    if "nScans" in selected.dimlabels:
        selected.mean("nScans")
    if selected.get_units(selected.dimlabels[-1]) != "s":
        selected.human_units()
    psd.plot(selected, "o", ax=ax4, alpha=0.5)
    ax4.set_title("table of integrals")
    psd.gridandtick(ax4)
    return selected, ax4
