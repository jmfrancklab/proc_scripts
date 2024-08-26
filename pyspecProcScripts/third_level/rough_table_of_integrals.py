from ..phasing import zeroth_order_ph, determine_sign, fid_from_echo
from ..simple_functions import select_pathway
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi


def rough_table_of_integrals(
    d, signal_range, fl=None, echo_like=False, title="", direct="t2"
):
    """manipulate s to generate a table of integrals (with only rough alignment)

    Returns
    =======
    s: nddata
        The table of integrals (collapse the direct dimension into a single number)
    ax3: Axes
        return the axis with the table of integrals, in case you want to add a fit!
    """
    assert fl is not None, "for now, fl can't be None"
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_figwidth(15)
    fig.set_figheight(6)
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle(title)
    # }}}
    # {{{ Apply overall zeroth order correction, select the pathway, and apply a rudimentary alignment
    d /= zeroth_order_ph(
        select_pathway(
            d[direct:signal_range].sum(direct), d.get_prop("coherence_pathway")
        )
    )
    if echo_like:
        # if it's echo like we will be applying fid_from_echo which requires
        # that the data passed still contains all coherence pathways
        s = select_pathway(
            d[direct:signal_range].C, d.get_prop("coherence_pathway")
        )
    else:
        s = select_pathway(
            d[direct:signal_range], d.get_prop("coherence_pathway")
        )
    # {{{ determine shift for rough alignment, but don't use this yet, because
    #     I want to see the signal resonance frequency
    s.ift(direct)
    shift = s * np.exp(-(pi**2) * s.fromaxis(direct) ** 2 * (2 * 50**2))
    shift.ft(direct)
    shift = shift[direct:signal_range].real.run(abs).argmax(direct)
    shift.set_error(None)
    # interestingly, if we were correcting here, we would multiply by the phase
    # here to perform the shift, which would be efficient!
    s.ft(direct)
    # }}}
    fl.image(
        s[direct:signal_range],
        ax=ax1,
        human_units=False,
    )
    ax1.set_title("Signal pathway / ph0")
    # }}}
    # {{{ Check phase variation along indirect
    mysign = determine_sign(
        s,
        signal_range,
    )
    s *= mysign
    fl.image(
        s[direct:signal_range],
        ax=ax2,
        human_units=False,
    )
    ax2.set_title("Check phase variation along indirect")
    # }}}
    if echo_like:
        d *= mysign
        s = fid_from_echo(d.set_error(None), d.get_prop("coherence_pathway"))
        s *= mysign
        s = select_pathway(
            s[direct:signal_range], d.get_prop("coherence_pathway")
        )
        fl.image(s[direct:signal_range], ax=ax3)
        ax3.set_title("FID sliced and phased")
    else:
        s *= mysign  # flip the sign back, so we get sensible integrals
    # {{{ generate the table of integrals
    # {{{ actually apply the rough alignment
    s.ift(direct)
    s *= np.exp(-1j * 2 * pi * shift * s.fromaxis(direct))
    s.ft(direct)
    # }}}
    s = s.real.integrate(direct).set_error(None)
    # }}}
    return s, ax3
