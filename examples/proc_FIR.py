"""
Process FIR experiment 
====================================================
Opens .h5 results file, uses rough_table_of_integrals() to roughly process dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
from Instruments.logobj import logobj
import numpy as np
import matplotlib.pyplot as plt

with psd.figlist_var() as fl:
    thisfile, exptype, post_proc, lookup = (
        "240924_13p5mM_TEMPOL_ODNP_1.h5",
        "ODNP_NMR_comp/ODNP",
        "spincore_IR_v3",
        prscr.lookup_table,
    )
    for nodename in [
        "FIR_noPower",
        "FIR_34dBm",
        ]:
        print(nodename)
        s = psd.find_file(
            thisfile,
            exp_type=exptype,
            expno=nodename,
            postproc=post_proc,
            lookup=prscr.lookup_table,
        )
        s.set_units("vd", "s")
        clock_correction = True
        indirect = "vd"
        direct = "t2"
        if clock_correction:
            # {{{ clock correction
            clock_corr = psd.nddata(np.linspace(-3, 3, 2500), "clock_corr")
            assert s.get_ft_prop(direct)
            if fl is not None:
                fl.next("before clock correction")
                fl.image(s)
            s_clock = prscr.select_pathway(s, s.get_prop("coherence_pathway")).mean('nScans').sum(
                direct
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
            plt.axvline(x=clock_corr, alpha=0.5, color="r")
            s *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
            s.ft(phase_dims)
            if fl is not None:
                fl.next("after auto-clock correction")
                fl.image(s)
            # }}}
        prscr.rough_table_of_integrals(s, fl=fl)
