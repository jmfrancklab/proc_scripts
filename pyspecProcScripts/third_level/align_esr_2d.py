import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from pyspecdata import pyspec_config

_figure_mode_setting = pyspec_config.get_setting(
    "figures", section="mode", environ="pyspecdata_figures"
)


def align_esr_2d(
    d,
    indirect_dim,
    direct_dim="$B_0$",
    correlation_slice=None,
    zerofill_factor=10,
    fl=None,
):
    r"""Correlation-align a 2D ESR data set in place.

    This is intended for ESR data where every spectrum shares a common direct
    field axis and the spectra are arranged along one indirect dimension, such
    as temperature or time.  Unlike :func:`align_esr`, this does not build a
    shared interpolation axis and does not apply any vertical scale matching.
    It only removes resonator-frequency drift by aligning each trace to the
    average spectrum, then subtracts the mean shift so the average field
    position is preserved.

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
    zerofill_factor: int default 10
        Zero-filling factor used only while calculating the correlation
        function.  This gives a finer shift axis without changing the number
        of points in the final aligned spectrum.
    fl: figlist_var default None
        If supplied, show diagnostic plots similar to :func:`align_esr`.

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
        fl.image(d.real)
        if _figure_mode_setting != "standard":
            fl.par_break()
        fl.next("2D ESR correlation", legend=True)
        if _figure_mode_setting != "standard":
            fl.par_break()
        fl.next("ESR shifts", legend=True)
        if _figure_mode_setting != "standard":
            fl.par_break()
        fl.next("Aligned 2D ESR", legend=True)

    baseline = d[direct_dim, -50:].C.mean(direct_dim)
    baseline_shape = []
    for this_dim in d.dimlabels:
        if this_dim in baseline.dimlabels:
            baseline_shape.append(d.data.shape[d.dimlabels.index(this_dim)])
        else:
            baseline_shape.append(1)
    d.data -= baseline.data.reshape(baseline_shape)
    d.set_ft_initial(direct_dim, "f")

    # Align each trace against the average spectrum.  The average gives a
    # stable series-wide reference without privileging any single indirect
    # point.
    average_spectrum = d.C.mean(indirect_dim)
    average_spectrum.ift(direct_dim)
    average_spectrum.run(np.conj)

    # Calculate all correlations at once.  The direct-axis Fourier transform
    # turns the product with the conjugated reference into a shift axis.
    # We zero fill only the correlation calculation, because this finds shifts
    # on a finer grid while leaving the final ESR spectra at their original
    # number of direct-axis points.
    correlation_pad = d.data.shape[d.dimlabels.index(direct_dim)]
    correlation_pad *= zerofill_factor
    correlation = d.C
    correlation.ift(direct_dim)
    correlation *= average_spectrum
    correlation.ft_new_startpoint(direct_dim, "f")
    correlation.ft(direct_dim, shift=True, pad=correlation_pad)
    if fl:
        fl.next("2D ESR correlation")
        fl.image(correlation.real)
    if correlation_slice is not None:
        correlation = correlation[direct_dim:correlation_slice]
    shifts = correlation.real.argmax(direct_dim)

    # Preserve the original average field position by forcing the mean shift to
    # zero, then store the full set of applied shifts for downstream review.
    shifts -= shifts.C.mean(indirect_dim).data
    d.set_prop("ESR_alignment_shifts", shifts)
    d.set_prop("ESR_alignment_direct_dim", direct_dim)
    d.set_prop("ESR_alignment_indirect_dim", indirect_dim)
    if fl:
        fl.next("ESR shifts")
        fl.plot(shifts, label="mean-preserving shift")
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
    if fl:
        fl.next("Aligned 2D ESR")
        fl.image(d)
    return d
