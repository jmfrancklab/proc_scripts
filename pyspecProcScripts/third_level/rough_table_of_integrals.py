from ..phasing import zeroth_order_ph, determine_sign, fid_from_echo
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
):
    """manipulate s to generate a table of integrals (with only rough alignment)

    Parameters
    ==========
    s : nddata
        Data with a single (dominant) peak, where you want to return the sign
        of the integral over all the data.
    signal_range : tuple (default None)
        Narrow slice range where signal resides.
        You probably want to get this from `find_peakrange` rather than trying
        to specify manually!
        (And expand the range that it gives you slightly)

        If this is set to None, it assumes a previous call to `find_peakrange` has set the
        `peakrange` property, and it uses that.
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
        Expand peakrange about its center by this much.

    Returns
    =======
    s : nddata
        The table of integrals (collapse the direct dimension into a single number).
        Processing is done in place.
    ax_last : Axes
        Return the axis with the table of integrals (plotted as `"o"`), in case you want to add a fit!
    """
    if signal_range is None:
        signal_range = s.getprop("peakrange")
    center_of_slice = np.mean(signal_range)
    signal_range_expanded = center_of_slice + expansion * r_[-0.5, 0.5] * np.diff(
        signal_range
    )
    assert fl is not None, "for now, fl can't be None"
    # {{{ set up subplots
    if echo_like:
        # If echo like we want an extra subplot to show the phased and FID sliced
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_figwidth(15)
    fig.set_figheight(6)
    fig.suptitle(title)
    fl.next("Raw Data with averaged scans", fig=fig)
    fl.skip_units_check()
    # }}}
    # To handle both older data where the coherence pathway was not set as a property
    # and handle newer data where it is set I have the following if/else
    signal_pathway = signal_pathway or s.get_prop("coherence_pathway")
    # {{{ Apply overall zeroth order correction, select the pathway, and apply a rudimentary alignment
    s /= zeroth_order_ph(
        select_pathway(s[direct:signal_range].sum(direct), signal_pathway)
    )
    s = select_pathway(s[direct:signal_range_expanded], signal_pathway)
    # {{{ determine shift for rough alignment
    s.ift(direct)
    shift = s * np.exp(-(pi**2) * s.fromaxis(direct) ** 2 * (2 * 50**2))
    shift.ft(direct)
    shift = shift.real.run(abs).argmax(direct)
    shift.set_error(None)
    # interestingly, if we were correcting here, we would multiply `s` by
    # the phase here to perform the shift, which would be efficient!
    s.ft(direct)
    # }}}
    fl.image(
        s,
        ax=ax1,
    )
    ax1.set_title("Signal pathway / ph0")
    # }}}
    # {{{ Check phase variation along indirect
    #     I don't want the data aligned to do this,  in case there is
    #     some type of second-order phase shift that becomes apparend
    mysign = determine_sign(
        s,
        signal_range,
    )
    s *= mysign
    fl.image(
        s,
        ax=ax2,
    )
    ax2.set_title("Check phase variation along indirect")
    # }}}
    if echo_like:
        signal_pathway = {}
        s = fid_from_echo(s.set_error(None), signal_pathway)
        s *= mysign
        fl.image(s, ax=ax3)
        ax3.set_title("FID sliced and phased")
    else:
        s *= mysign  # flip the sign back, so we get sensible integrals
    # {{{ generate the table of integrals
    # {{{ actually apply the rough alignment
    #     Note that we center about the center of the slice, so our
    #     sliced integration works out OK, not about zero.
    s.ift(direct)
    s *= np.exp(-1j * 2 * pi * (shift - center_of_slice) * s.fromaxis(direct))
    s.ft(direct)
    # }}}
    s = s[direct:signal_range].real.integrate(direct).set_error(None)
    ax_last = ax4 if echo_like else ax3
    fl.plot(s, "o", ax=ax_last)
    psd.gridandtick(plt.gca())
    ax_last.set_title("table of integrals")
    # }}}
    return s, ax_last
