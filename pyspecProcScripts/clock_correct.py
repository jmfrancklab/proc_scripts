import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt
import numpy as np


def clock_correct(s, fl=None, indirect="vd"):
    """Do clock correction.

    On spincore, the Tx carrier (clock) and Rx carrier do not seem to be
    exactly the same.
    This can be proven with a simple experiment that involves changing
    the time between phase reset and the start of the pulses.
    It is subtle, so it only affects things with longer delays, like T₁
    measurements.

    Parameters
    ==========
    s : nddata
        The data to be corrected.

    indirect : str
        The dimension along which to clock correct.

    fl : figlist

    Returns
    ==========
    s: nddata
        modified in-place

    """
    clock_corr = psd.nddata(np.linspace(-3, 3, 2500), "clock_corr")
    assert s.get_ft_prop("t2")
    if fl is not None:
        fl.next("before clock correction")
        fl.image(s)
    s_clock = (
        prscr.select_pathway(s, s.get_prop("coherence_pathway"))
        .mean("nScans")
        .sum("t2")
    )
    phase_dims = [j for j in s.dimlabels if j.startswith("ph")]
    s.ift(phase_dims)
    min_index = abs(s_clock).argmin(indirect, raw_index=True).item()
    s_clock *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
    s_clock[indirect, : min_index + 1] *= -1
    s_clock.sum(indirect).run(abs)
    if fl is not None:
        fl.next("clock correction")
        fl.plot(s_clock, ".", alpha=0.7)
    clock_corr = s_clock.argmax("clock_corr").item()
    if fl is not None:
        plt.axvline(x=clock_corr, alpha=0.5, color="r")
    s *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
    s.ft(phase_dims)
    if fl is not None:
        fl.next("after auto-clock correction")
        fl.image(s)
    return s
