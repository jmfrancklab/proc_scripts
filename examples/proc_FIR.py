"""
Process FIR experiment 
====================================================
Opens .h5 results file, uses :func:`rough_table_of_integrals` to roughly process
dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
import numpy as np
import sympy
import matplotlib.pyplot as plt
from itertools import cycle

T1_list = []
clock_correction = True

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
        fl.basename = nodename # this is a good example of how we can
        #                        use basename to make it easy to deal
        #                        with multiple datasets
        s = psd.find_file(
            thisfile,
            exp_type=exptype,
            expno=nodename,
            postproc=post_proc,
            lookup=prscr.lookup_table,
        )
        indirect = "vd"
        direct = "t2"
        if clock_correction:
            # {{{ clock correction
            s = prscr.clock_correct(s, fl=fl)
            # }}}
        s, dropped_dims = s.squeeze(return_dropped=True)
        s, ax_last = prscr.rough_table_of_integrals(s, fl=fl)
        # TODO ☐: I'm assuming this is copied from previous.  I think that this
        # is probably something we want to address in rough_table_of_integrals
        # (if we have multiple scans, we want the result to be signal
        # averaged).  Likely you would signal average after alignment, but
        # before you actually take the integral over range of frequencies.
        # But, I think it's useful for you to think about this as well and
        # implement a solution.
        if "nScans" in s.dimlabels:
            s.mean("nScans")
        # TODO ☐: below, you are tearing apart the nddata object and putting it
        # back together -- that's always a bad sign
        Mi, R1, vd = sympy.symbols("M_inf R_1 vd", real=True)
        psd.logger.debug(psd.strm("acq keys", s.get_prop("acq_params")))
        W = (
            s.get_prop("acq_params")["FIR_rep"] * 1e-6
            + s.get_prop("acq_params")["acq_time_ms"] * 1e-3
        )
        # (delete when done) the code didn't match our definition of "clean"
        # because you're creating many variable names for the same dataset.
        # Also, you are defining variables in order to use them only on one
        # other line of code, rather than combining.  The way around a lot of
        # this is to think about how you can *reorder* the code.
        s = psd.lmfitdata(s)
        s.functional_form = Mi * (
            1 - (2 - sympy.exp(-W * R1)) * sympy.exp(-vd * R1)
        )
        # TODO ☐: the M_inf on the next line should be set based on the data
        # (some multiple of the max and min)
        s.set_guess(
            M_inf=dict(value=3.9e4, min=0, max=2e7),
            R_1=dict(value=0.8, min=0.01, max=100),
        )
        s.fit()
        s_fit = s.eval(200)
        psd.plot(s_fit,ax=ax_last, alpha=0.5) # here, we plot the fit
        #                                       together with the table of
        #                                       integrals.
        # TODO ☐: next pull from code below to collect the normalized fits into a single plot.
        # I don't really see the value of the unnormalized anymore because we have each
        # individual curve together with the table of integrals, now.
        if False: # JF has not reviewed this -- needs to be re-written
            #       consistently w/ above.  Stuff that's not used can
            #       just be removed
            fl.next("IR fit - before norm")
            fl.plot(s_int, "o", label="%s" % nodename)
            f.fit()
            T1 = 1.0 / f.output("R_1")
            Mi = f.output("M_inf")
            T1_list.append(T1)
            f.set_units("vd", "s")
            fit = f.eval(100)
            fit.set_plot_color(f.get_plot_color())
            fl.plot(
                fit,
                ls="-",
                alpha=0.5,
                label="fit for %s" % nodename,
            )
            fl.next("IR fit - before norm - %s" % nodename)
            fl.plot(s_int, "o", label="%s" % nodename)
            fl.plot(
                fit,
                ls="-",
                alpha=0.5,
                label="fit for %s" % nodename,
            )
            ax = plt.gca()
