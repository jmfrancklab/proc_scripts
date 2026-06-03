import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from pyspecdata import nddata, pyspec_config

from ..cumulant_rms_func import cumulant_rms

_figure_mode_setting = pyspec_config.get_setting(
    "figures", section="mode", environ="pyspecdata_figures"
)


def calculate_esr_2d_shift(
    d,
    indirect_dim,
    direct_dim="$B_0$",
    correlation_slice=None,
    zerofill_factor=200,
    ax=None,
):
    """Calculate uncentered ESR shifts for one locally similar chunk.

    The caller is responsible for selecting a chunk where spectra have similar
    lineshapes.  This routine calculates the correlation against the local
    average spectrum, optionally plots that correlation on an existing Axes,
    and returns uncentered shifts from the zero-filled correlation axis.
    """
    # The local average is used as the reference so that each chunk compares
    # spectra only against nearby, similar lineshapes.
    average_spectrum = d.C.mean(indirect_dim)
    average_spectrum.ift(direct_dim)
    average_spectrum.run(np.conj)

    # Zero fill only the correlation calculation.  This refines the shift grid
    # without changing the number of points in the aligned ESR spectra.
    correlation_pad = d.data.shape[d.dimlabels.index(direct_dim)]
    correlation_pad *= zerofill_factor
    correlation = d.C
    correlation.ift(direct_dim)
    correlation *= average_spectrum
    correlation.ft_new_startpoint(direct_dim, "f")
    correlation.ft(direct_dim, shift=True, pad=correlation_pad)
    if ax is not None:
        correlation_for_plot = correlation.real.C.reorder(
            [indirect_dim, direct_dim]
        )
        ax.imshow(
            correlation_for_plot.data,
            aspect="auto",
            origin="lower",
        )
        ax.set_title("correlation")
        ax.set_xlabel(direct_dim)
        ax.set_ylabel(indirect_dim)
    if correlation_slice is not None:
        correlation = correlation[direct_dim:correlation_slice]
    return correlation.real.argmax(direct_dim)


def align_esr_2d(
    d,
    indirect_dim,
    direct_dim="$B_0$",
    correlation_slice=None,
    zerofill_factor=200,
    n_chunks=5,
    chunk_extender_frac=0.1,
    fl=None,
):
    r"""Correlation-align a 2D ESR data set in place.

    This is intended for ESR data where every spectrum shares a common direct
    field axis and the spectra are arranged along one indirect dimension, such
    as temperature or time.  Unlike :func:`align_esr`, this does not build a
    shared interpolation axis and does not apply any vertical scale matching.
    It only removes resonator-frequency drift by aligning each trace to local
    average spectra from regions with similar cumulant RMS values, then
    subtracts the mean shift so the average field position is preserved.

    Parameters
    ==========
    d: nddata
        The 2D ESR data set to adjust in place.
    indirect_dim: str
        The dimension indexing the ESR spectra to be aligned.
    direct_dim: str, default ``"$B_0$"``
        The ESR field dimension.
    correlation_slice: tuple default None
        Optional range of direct-axis shifts over which to search for the
        correlation maximum.
    zerofill_factor: int default 200
        Zero-filling factor used only while calculating the correlation
        function.  This gives a finer shift axis without changing the number
        of points in the final aligned spectrum.
    n_chunks: int default 5
        Number of cumulant-RMS chunks used for local correlation references.
    chunk_extender_frac: float default 0.1
        Fraction of a chunk width used to extend both chunk ends.  The default
        gives 20% overlap between neighboring chunks.
    fl: figlist_var default None
        If supplied, show diagnostic plots for the raw data, correlation
        chunks, and stitched ESR shift curve.

    Returns
    =======
    nddata
        The same object supplied as ``d``, baseline corrected, shifted in
        place, and converted to its real component.
    """
    # Remove a simple right-edge baseline from each spectrum before alignment,
    # matching the lightweight baseline correction used by align_esr.
    if fl:
        if _figure_mode_setting != "standard":
            fl.par_break()
        fl.next("Raw 2D ESR", legend=True)
        fl.image(d.real.C.reorder([indirect_dim, direct_dim]))
        if _figure_mode_setting != "standard":
            fl.par_break()

    baseline = d[direct_dim, -50:].C.mean(direct_dim)
    baseline_shape = []
    for this_dim in d.dimlabels:
        if this_dim in baseline.dimlabels:
            baseline_shape.append(d.data.shape[d.dimlabels.index(this_dim)])
        else:
            baseline_shape.append(1)
    d.data -= baseline.data.reshape(baseline_shape)
    d.set_ft_initial(direct_dim, "f")

    # Use cumulant RMS to choose chunks that span similar amounts of spectral
    # change, rather than simply equal numbers of indirect points.
    indirect_axis = d.getaxis(indirect_dim)
    n_indirect = d.data.shape[d.dimlabels.index(indirect_dim)]
    cumulant = cumulant_rms(d, indirect_dim, direct_dim)
    cumulant_data = cumulant.data.real.ravel()
    if cumulant_data[-1] == 0:
        cumulant_data = np.linspace(0, 1, len(cumulant_data))
    else:
        cumulant_data /= cumulant_data[-1]

    if fl:
        n_rows = int(np.ceil(np.sqrt(n_chunks)))
        n_cols = int(np.ceil(float(n_chunks) / n_rows))
        fig, ax_array = plt.subplots(n_rows, n_cols)
        fl.next("2D ESR correlation chunks", fig=fig)
        ax_array = np.asarray(ax_array).ravel()
        if _figure_mode_setting != "standard":
            fl.par_break()
    else:
        ax_array = [None] * n_chunks

    chunk_shift_plots = []
    shift_sum = np.zeros(n_indirect)
    shift_count = np.zeros(n_indirect)
    for j in range(n_chunks):
        # The extender is divided by n_chunks, so the absolute overlap becomes
        # smaller when fewer chunks are requested.
        start_frac = float(j - chunk_extender_frac) / n_chunks
        stop_frac = float(j + 1 + chunk_extender_frac) / n_chunks
        start_frac = max(0, start_frac)
        stop_frac = min(1, stop_frac)
        if cumulant.data.size > 1:
            start_idx = np.searchsorted(cumulant_data, start_frac, side="left")
            stop_idx = np.searchsorted(cumulant_data, stop_frac, side="right")
            stop_idx += 1
        else:
            start_idx = int(np.floor(start_frac * n_indirect))
            stop_idx = int(np.ceil(stop_frac * n_indirect))
        start_idx = max(0, min(start_idx, n_indirect - 1))
        stop_idx = max(start_idx + 1, min(stop_idx, n_indirect))

        this_chunk = d[indirect_dim, start_idx:stop_idx]
        this_shift = calculate_esr_2d_shift(
            this_chunk,
            indirect_dim,
            direct_dim=direct_dim,
            correlation_slice=correlation_slice,
            zerofill_factor=zerofill_factor,
            ax=ax_array[j],
        )
        if ax_array[j] is not None:
            ax_array[j].set_title(
                "%s %0.3g to %0.3g"
                % (
                    indirect_dim,
                    indirect_axis[start_idx],
                    indirect_axis[stop_idx - 1],
                )
            )

        # Align the uncentered chunk shifts to the already-stitched curve using
        # their overlap, then accumulate and average all overlapped points.
        overlap_mask = shift_count[start_idx:stop_idx] > 0
        if np.any(overlap_mask):
            stitched_so_far = shift_sum[start_idx:stop_idx][overlap_mask]
            stitched_so_far /= shift_count[start_idx:stop_idx][overlap_mask]
            this_shift.data -= np.mean(
                this_shift.data[overlap_mask] - stitched_so_far
            )
        shift_sum[start_idx:stop_idx] += this_shift.data
        shift_count[start_idx:stop_idx] += 1
        chunk_shift_plots.append(this_shift.C)

    if fl:
        for j in range(n_chunks, len(ax_array)):
            ax_array[j].axis("off")
        fig.tight_layout()

    full_shift_data = shift_sum / shift_count
    center_shift = np.mean(full_shift_data)
    full_shift_data -= center_shift
    for this_shift in chunk_shift_plots:
        this_shift.data -= center_shift

    shifts = nddata(full_shift_data, [n_indirect], [indirect_dim])
    shifts.setaxis(indirect_dim, indirect_axis)
    d.set_prop("ESR_alignment_shifts", shifts)
    d.set_prop("ESR_alignment_direct_dim", direct_dim)
    d.set_prop("ESR_alignment_indirect_dim", indirect_dim)
    d.set_prop("ESR_alignment_n_chunks", n_chunks)
    d.set_prop("ESR_alignment_zerofill_factor", zerofill_factor)
    if fl:
        fl.next("ESR shifts", legend=True)
        for j, this_shift in enumerate(chunk_shift_plots):
            plt.plot(
                this_shift.getaxis(indirect_dim),
                this_shift.data,
                alpha=0.5,
                label="chunk %d" % (j + 1),
            )
        plt.plot(
            shifts.getaxis(indirect_dim),
            shifts.data,
            "x",
            alpha=0.5,
            color="k",
            label="stitched average",
        )
        plt.ylabel(direct_dim)

    # Apply the shifts in place using the same Fourier phase-ramp convention
    # as align_esr, then return to the direct domain and keep only real data.
    d.ift(direct_dim)
    direct_axis_shape = []
    shift_shape = []
    for this_dim in d.dimlabels:
        if this_dim == direct_dim:
            direct_axis_shape.append(d.data.shape[d.dimlabels.index(this_dim)])
            shift_shape.append(1)
        elif this_dim in shifts.dimlabels:
            direct_axis_shape.append(1)
            shift_shape.append(d.data.shape[d.dimlabels.index(this_dim)])
        else:
            direct_axis_shape.append(1)
            shift_shape.append(1)
    d.data *= np.exp(
        -1j
        * 2
        * pi
        * d.fromaxis(direct_dim).data.reshape(direct_axis_shape)
        * shifts.data.reshape(shift_shape)
    )
    d.ft(direct_dim)
    d.run(np.real)
    return d
