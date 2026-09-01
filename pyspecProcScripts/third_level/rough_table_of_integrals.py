from ..phasing import (
    zeroth_order_ph,
    determine_sign,
    fid_from_echo,
    det_inh_bounds,
)
from ..simple_functions import select_pathway
import matplotlib.pyplot as plt
import numpy as np
from numpy import r_
from numpy import pi
import pyspecdata as psd


def rough_table_of_integrals(
    s,
    signal_range=None,
    signal_pathway=None,
    fl=None,
    echo_like=True,
    title="",
    direct="t2",
    expansion=2,
    peak_lowest_thresh=0.1,
    inc_plot_color=True,
):
    """Manipulate ``s`` to generate a table of integrals
    (with only rough alignment)

    Parameters
    ==========
    s : nddata
        Data with a single (dominant) peak, where you want to return the sign
        of the integral over all the data.
    signal_range : tuple (default None)
        Narrow integration bounds where the absorptive signal resides.
        You probably want to get this from `det_inh_bounds` rather than trying
        to specify manually!

        If this is set to the string `'inh_bounds'`,
        it assumes a previous call to `det_inh_bounds`
        has set the
        `'inh_bounds'` property, and it uses that.

        If this is set to None, then it assumes you want to call
        `det_inh_bounds`
    signal_pathway : dict (default None)
        If None, the function will go into the properties of the data looking
        for the "coherence_pathway" property.
        If that doesn't exist the user needs to feed a dictionary containing
        the coherence transfer pathway where the signal resides.
        If the signal has already been sliced along the coherence dimensions,
        you can feed an empty dictionary (note the *very* different behavior
        *vs* `None`)
    echo_like : boolean
        If the signal has an echo the `fid_from_echo` function will be applied
        prior to integration to apply the FID slice as well as proper phasing
        and peak selection.
    title : str
        Title for the returned figure.
    direct : str (default "t2")
        Name of direct dimension.
    expansion : float (default 2)
        Expand ``inh_bounds`` about its center by this much only when locating
        the rough frequency shifts.  This does not set the echo-processing or
        integration bounds.
    peak_lowest_thresh : float
        Lowest threshold passed to :func:`det_inh_bounds`.
    inc_plot_color: boolean
        assume that we are processing multiple datasets, and want to increment
        the color counter with every run of this function.

    Returns
    =======
    s : nddata
        The table of integrals (collapse the direct dimension into a single
        number).
        Processing is done in place.
    ax4 : Axes or None
        Return the axis with the table of integrals (plotted as `"o"`), in case
        you want to add a fit.  Returns None when ``fl`` is None.
    """
    if signal_pathway is None:
        signal_pathway = s.get_prop("coherence_pathway")
    if signal_pathway is None:
        raise ValueError(
            "signal_pathway must be supplied or stored as a data property"
        )
    if signal_range is None:
        det_inh_bounds(
            s,
            peak_lowest_thresh,
            fl=fl,
            direct=direct,
            echo_like=echo_like,
            signal_pathway=signal_pathway,
        )
        signal_range = s.get_prop("inh_bounds")
    elif isinstance(signal_range, str):
        if signal_range == "inh_bounds":
            signal_range = s.get_prop("inh_bounds")
        else:
            raise ValueError(
                "signal_range strings must contain only 'inh_bounds'"
            )
    signal_range = np.asarray(signal_range, dtype=float)
    if signal_range.shape != (2,) or signal_range[0] >= signal_range[1]:
        raise ValueError("signal_range must contain two ascending bounds")
    center_of_range = np.mean(signal_range)
    half_range = 0.5 * np.diff(signal_range).item()
    signal_range_expanded = (
        center_of_range + expansion * r_[-1, 1] * half_range
    )
    ax4 = None
    # {{{ set up subplots
    if fl is not None:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        fl.next("Raw Data with averaged scans", fig=fig)
        fig.suptitle(title)
        fl.skip_units_check()
    # }}}
    # {{{ Apply overall zeroth order correction, select the pathway, and apply
    #     a rudimentary alignment
    s /= zeroth_order_ph(
        select_pathway(s[direct:signal_range].sum(direct), signal_pathway)
    )
    s = select_pathway(s, signal_pathway)
    # {{{ determine shift for rough alignment
    # Locate the shifts on a temporary narrow copy, while preserving the full
    # spectrum for linewidth-aware echo processing below.
    shift = s[direct:signal_range_expanded].C
    shift.ift(direct)
    shift *= np.exp(-(pi**2) * shift.fromaxis(direct) ** 2 * (2 * 50**2))
    shift.ft(direct)
    shift = shift.real.run(abs).argmax(direct)
    shift.set_error(None)
    # interestingly, if we were correcting here, we would multiply `s` by
    # the phase here to perform the shift, which would be efficient!
    # }}}
    if fl is not None:
        fl.image(
            s,
            ax=ax1,
            interpolation="auto",
        )
        ax1.set_title("extract signal pathway")
    # }}}
    # {{{ Check phase variation along indirect
    #     I don't want the data aligned to do this,  in case there is
    #     some type of second-order phase shift that becomes apparend
    mysign = determine_sign(
        s,
        signal_range,
    )
    s *= mysign
    if fl is not None:
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
            inh_bounds=signal_range,
            max_alignment_shift=abs(shift - center_of_range).data.max().item(),
            fl=fl,
        )
        s *= mysign
    else:
        s *= mysign  # flip the sign back, so we get sensible integrals
    # {{{ generate the table of integrals
    # {{{ actually apply the rough alignment
    #     Note that we center about the center of the slice, so our
    #     sliced integration works out OK, not about zero.
    s.ift(direct)
    s *= np.exp(-1j * 2 * pi * (shift - center_of_range) * s.fromaxis(direct))
    s.ft(direct)
    if fl is not None:
        fl.image(
            s,
            ax=ax3,
            interpolation="auto",
        )
        ax3.set_title(
            "FID sliced" + (", phased," if echo_like else "") + " and aligned"
        )
    # }}}
    s = s[direct:signal_range].real.integrate(direct).set_error(None)
    if inc_plot_color:
        s.set_plot_color_next()
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    if s.get_units(s.dimlabels[-1]) != "s":
        s.human_units()  # because the units aren't under the control of
        #                  the figure list, we're going to convert to
        #                  "human" units here
    if fl is not None:
        psd.plot(s, "o", ax=ax4, alpha=0.5)
        ax4.set_title("table of integrals")
        psd.gridandtick(ax4)
    # }}}
    return s, ax4
