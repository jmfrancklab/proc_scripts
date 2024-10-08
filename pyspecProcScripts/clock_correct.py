import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt
import numpy as np

def clock_correct(s, fl=None):
    """ Do clock correction
    
    Parameters
    ==========
    s:

    fl: figlist_var()
    
    Returns
    ==========
    s:

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
    min_index = abs(s_clock).argmin("vd", raw_index=True).item()
    s_clock *= np.exp(-1j * clock_corr * s.fromaxis("vd"))
    s_clock["vd", : min_index + 1] *= -1
    s_clock.sum("vd").run(abs)
    if fl is not None:
        fl.next("clock correction")
        fl.plot(s_clock, ".", alpha=0.7)
    clock_corr = s_clock.argmax("clock_corr").item()
    plt.axvline(x=clock_corr, alpha=0.5, color="r")
    s *= np.exp(-1j * clock_corr * s.fromaxis("vd"))
    s.ft(phase_dims)
    if fl is not None:
        fl.next("after auto-clock correction")
        fl.image(s)
    return s
