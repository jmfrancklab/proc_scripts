from ..correlation_alignment import correl_align
from ..calc_error import calc_masked_variance
from ..envelope import fit_envelope, L2G
from ..phasing import (
    zeroth_order_ph,
    determine_sign,
    fid_from_echo,
    find_peakrange,
    hermitian_function_test,
)
from ..simple_functions import select_pathway
import matplotlib.pyplot as plt
import numpy as np
from numpy import r_
from numpy import pi
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
    max_shift=300.0,
    center_aligned_peak=True,
    show_alignment_diagnostics=False,
    propagate_error=True,
    excluded_pathways=None,
    require_unitary_error=True,
    equal_energy_apodization=True,
    apodization_lambda=None,
    show_apodization_diagnostics=False,
):
    """Generate a table of integrals after correlation alignment.

    This follows the same processing path as ``rough_table_of_integrals``:
    find/select the signal range, apply zeroth-order phase correction, check
    phase variation, optionally slice the FID from an echo, integrate the real
    signal, average scans, and plot the resulting table.  The difference is
    that the frequency shift is determined by ``correl_align`` instead of the
    rough Gaussian-filtered argmax alignment.

    Parameters
    ==========
    s : nddata
        Data with a single dominant peak.  The direct dimension is expected to
        start in the frequency domain, and phase-cycling dimensions should be
        in the coherence-transfer domain.
    signal_range : tuple (default None)
        Narrow slice range where signal resides.  If None, call
        ``find_peakrange``.  If set to ``"peakrange"``, use the existing
        ``peakrange`` property.
    signal_pathway : dict (default None)
        Coherence transfer pathway for the signal.  If None, read the
        ``coherence_pathway`` property.
    fl : figlist
        Required figure list.
    echo_like : boolean
        If true, call ``fid_from_echo`` before integration.
    title : str
        Figure title.
    direct : str
        Direct dimension.
    expansion : float
        Expand the peak range about its center before plotting/checking phase.
    peak_lower_thresh : float
        Passed to ``find_peakrange``.
    inc_plot_color : boolean
        If true, increment the plot color after generating the table.
    repeat_dims : str or list (default None)
        Dimensions that contain repeated signal traces and should be aligned.
        If None, infer all non-phase, non-direct dimensions except ``nScans``
        and ``non_repeat_dims``.
    non_repeat_dims : str or list (default None)
        Dimensions that should not be treated as repeated traces.
    alignment_sigma : float
        Width parameter for the frequency-domain Gaussian mask used by
        ``correl_align``.
    max_shift : float
        Maximum frequency shift allowed by ``correl_align``.
    center_aligned_peak : boolean
        If true, recenter the integration window on the peak of the aligned
        average.  This is useful for FIR nodes, where the aligned average is a
        good representative spectrum.  Set false for E(p) power trajectories,
        where a single averaged center can be dominated by high-SNR/high-power
        points.
    show_alignment_diagnostics : boolean
        If true, send the figure list into the Hermitian centering and
        correlation-alignment diagnostics.  The main four-panel table plot is
        always generated through ``fl``.
    propagate_error : boolean
        If true, estimate spectral-point variance from non-selected coherence
        pathways and propagate it through the direct-dimension integral.
    excluded_pathways : list of dicts
        Pathways to mask before estimating the variance.  By default, only the
        selected signal pathway is masked, so signal/artifact leakage into
        other pathways contributes to the propagated error.
    require_unitary_error : bool
        Passed to ``calc_masked_variance``.  Keep true for freshly processed
        data with unitary phase-cycling Fourier transforms; set false for
        legacy/postprocessed data where the coherence-domain variance is used
        with its existing scaling.
    equal_energy_apodization : bool
        If true, apply the DCCT equal-energy Lorentzian-to-Gaussian
        apodization in the time domain after Hermitian centering and before
        correlation alignment.
    apodization_lambda : float or None
        Lorentzian linewidth used for the L2G apodization.  If None, estimate
        it from the selected pathway envelope.
    show_apodization_diagnostics : bool
        If true, show envelope-fit diagnostic plots used to determine the
        apodization linewidth.

    Returns
    =======
    s : nddata
        The table of integrals.
    ax4 : Axes
        Axis containing the table plot.
    """
    if signal_range is None:
        center_of_range, half_range = find_peakrange(
            s, fl=fl, direct=direct, peak_lower_thresh=peak_lower_thresh
        )
        signal_range = s.get_prop("peakrange")
    else:
        if signal_range == "peakrange":
            signal_range = s.get_prop("peakrange")
        center_of_range = np.mean(signal_range)
        half_range = 0.5 * np.diff(signal_range).item()
    signal_range_expanded = (
        center_of_range + expansion * r_[-1, 1] * half_range
    )
    assert fl is not None, "for now, fl can't be None"
    signal_pathway = (
        s.get_prop("coherence_pathway")
        if signal_pathway is None
        else signal_pathway
    )
    if signal_pathway is None:
        raise ValueError(
            "table_of_integrals needs signal_pathway or a coherence_pathway"
            " property so correlation alignment can mask pathways"
        )
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
            " correlation alignment; pass repeat_dims explicitly if inference"
            " is not correct"
        )

    # {{{ set up subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.suptitle(title)
    fl.next("Raw Data with correlation alignment", fig=fig)
    fl.skip_units_check()
    # }}}

    # {{{ Apply zeroth-order phase and correlation alignment
    s /= zeroth_order_ph(
        select_pathway(s[direct:signal_range].sum(direct), signal_pathway)
    )
    s.ift(direct)
    s[direct] -= s.getaxis(direct)[0]
    s_for_center = s.C
    for thisdim in repeat_dims + ["nScans"]:
        if thisdim in s_for_center.dimlabels:
            s_for_center.mean(thisdim)
    best_shift = hermitian_function_test(
        select_pathway(s_for_center, signal_pathway),
        direct=direct,
        fl=fl if show_alignment_diagnostics else None,
    )
    s.setaxis(direct, lambda x: x - best_shift).register_axis({direct: 0})
    if equal_energy_apodization:
        # {{{ Apply DCCT equal-energy apodization after Hermitian centering
        # L2G is a two-sided time-domain transformation, so the origin matters.
        # Apply it only after the Hermitian test has located the echo/FID
        # origin; otherwise the window is centered on the acquisition start and
        # can create two apparent frequency centers in E(p).
        if apodization_lambda is None:
            apodization_lambda = fit_envelope(
                select_pathway(s.C, signal_pathway),
                direct=direct,
                fl=fl if show_apodization_diagnostics else None,
            )
        s *= L2G(apodization_lambda, criterion="energy")(s.fromaxis(direct))
        # }}}
    s.ft(direct)
    mysign_for_alignment = (
        select_pathway(s, signal_pathway).C.real.sum(direct).run(np.sign)
    )
    mysign_for_alignment[lambda x: x == 0] = 1

    def frq_mask(x):
        nu_center = (
            select_pathway(x, signal_pathway).mean("repeats").argmax(direct)
        )
        return x * np.exp(
            -((x.fromaxis(direct) - nu_center) ** 2)
            / (4 * alignment_sigma**2)
        )

    def coherence_unmask(coh_array):
        if len(signal_pathway) == 0:
            coh_array.data[:] = 1
            return coh_array
        for pathway_dict in (signal_pathway, {k: 0 for k in signal_pathway}):
            thisslice = coh_array
            for j, (k, v) in enumerate(pathway_dict.items()):
                if j == len(pathway_dict) - 1:
                    thisslice[k, v] = 1
                else:
                    thisslice = thisslice[k, v]
        return coh_array

    # {{{ Align repeated measurements to their optimized common average
    opt_shift = correl_align(
        s * mysign_for_alignment,
        frq_mask_fn=frq_mask,
        coherence_unmask_fn=coherence_unmask,
        repeat_dims=repeat_dims,
        non_repeat_dims=non_repeat_dims,
        max_shift=max_shift,
        direct=direct,
        fl=fl if show_alignment_diagnostics else None,
    )
    # }}}
    s.ift(direct).ift(list(signal_pathway.keys()))
    s *= np.exp(-1j * 2 * pi * opt_shift * s.fromaxis(direct))
    s.ft(direct).ft(list(signal_pathway.keys()))
    aligned_for_error = s.C
    if center_aligned_peak:
        # {{{ Recenter the integration window on the aligned peak
        # The correlation step aligns repeated traces to each other; it does
        # not force the absolute resonance position to zero.  After the
        # relative alignment, find the peak in the aligned average and rebuild
        # the signal window around that frequency so the plotted/integrated
        # slice is centered on the aligned signal rather than on the
        # pre-alignment estimate.
        aligned_peak = abs(
            select_pathway(s[direct:signal_range_expanded].C, signal_pathway)
        )
        aligned_peak.mean_all_but([direct])
        center_of_range = aligned_peak.argmax(direct).item()
        signal_range = center_of_range + r_[-1, 1] * half_range
        signal_range_expanded = (
            center_of_range + expansion * r_[-1, 1] * half_range
        )
        # }}}
    s = select_pathway(s[direct:signal_range_expanded], signal_pathway)
    fl.image(
        s,
        ax=ax1,
        interpolation="auto",
    )
    ax1.set_title("extract signal pathway\nand correlation align")
    # }}}

    # {{{ Check phase variation along indirect
    mysign = determine_sign(
        s,
        signal_range,
    )
    s *= mysign
    integral_error = None
    spectral_datapoint_variance = None
    if propagate_error and len(signal_pathway) > 0:
        # {{{ Propagate DCCT leakage/noise from non-selected pathways
        # The DCCT paper motivates inspecting all coherence-transfer pathways
        # rather than treating them as discarded data.  Here, the selected
        # signal pathway is masked, while the same frequency window in all
        # other pathways is retained.  Therefore signal leakage and artifacts
        # that land in non-selected pathways increase the spectral-point
        # variance that is propagated into the final integral error.
        error_source = aligned_for_error[direct:signal_range_expanded]
        error_source *= mysign
        if echo_like:
            error_source = fid_from_echo(
                error_source.set_error(None),
                signal_pathway,
                frq_center=center_of_range,
                frq_half=half_range,
            )
        error_source = error_source[direct:signal_range]
        error_source.set_prop("coherence_pathway", signal_pathway)
        error_indirect_dims = [
            j
            for j in error_source.dimlabels
            if j != direct and not j.startswith("ph")
        ]
        spectral_datapoint_variance = calc_masked_variance(
            error_source,
            excluded_frqs=[],
            indirect=error_indirect_dims,
            excluded_pathways=excluded_pathways,
            direct=direct,
            require_unitary=require_unitary_error,
        )
        # }}}
    fl.image(
        s,
        ax=ax2,
        interpolation="auto",
    )
    ax2.set_title("check phase variation\nalong indirect")
    # }}}

    if echo_like:
        signal_pathway = {}
        s = fid_from_echo(
            s.set_error(None),
            signal_pathway,
            frq_center=center_of_range,
            frq_half=half_range,
            fl=fl if show_alignment_diagnostics else None,
        )
        s *= mysign
    else:
        s *= mysign

    # {{{ generate the table of integrals
    fl.image(
        s,
        ax=ax3,
        interpolation="auto",
    )
    ax3.set_title(
        "FID sliced" + (", phased," if echo_like else "") + " and aligned"
    )
    final_slice = s[direct:signal_range]
    if final_slice.shape[direct] < 2:
        # {{{ Recover from an empty final slice after echo/FID conversion
        # The alignment window is chosen in the original frequency-domain
        # data.  For echo-like data, fid_from_echo can rebuild the direct axis
        # so the old absolute frequency window no longer contains enough
        # points.  Keep the same width, but center the final integration on
        # the peak that is present in the aligned/FID-converted data.
        final_peak = abs(s.C)
        final_peak.mean_all_but([direct])
        center_of_range = final_peak.argmax(direct).item()
        signal_range = center_of_range + r_[-1, 1] * half_range
        final_slice = s[direct:signal_range]
        if final_slice.shape[direct] < 2:
            raise ValueError(
                "table_of_integrals could not find enough points in the "
                "final integration window"
            )
        # }}}
    s = final_slice.real
    if spectral_datapoint_variance is not None:
        df = s.get_ft_prop(direct, "df")
        if df is None:
            df = np.diff(s.getaxis(direct)[r_[0, 1]]).item()
        N = s.shape[direct]
        integral_error = (spectral_datapoint_variance.C * df**2 * N).run(
            np.sqrt
        )
    s.integrate(direct).set_error(None)
    if integral_error is not None:
        if set(integral_error.dimlabels) == set(s.dimlabels):
            integral_error.reorder(s.dimlabels)
        s.set_error(integral_error.data)
    if inc_plot_color:
        s.set_plot_color_next()
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    if s.get_units(s.dimlabels[-1]) != "s":
        s.human_units()
    psd.plot(s, "o", ax=ax4, alpha=0.5)
    ax4.set_title("table of integrals")
    psd.gridandtick(ax4)
    # }}}
    return s, ax4
