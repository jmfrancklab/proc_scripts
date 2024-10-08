"""
Process FIR experiment 
====================================================
Opens .h5 results file, uses :func:`rough_table_of_integrals` to roughly
process
dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy
import matplotlib.pyplot as plt

T1_list = []
clock_correction = True
plot_fit = True
# TODO ☐: to directly answer your question, the thing you get from
#         matplotlib isn't a cycle, but a "cycler", which is different.
#         You need to feed it to "cycle" from "itertools" in order to
#         make it into a "cycle" and be able to use "next" on it.
#         BUT you should't do that, anyways.  The nddata object that
#         comes out of rough_table_of_integrals will always plot with
#         the same color (as will copies of it).  You can also get the
#         color of that object (e.g. for use in an unrelated dataset
#         axvline, etc) with `s.get_plot_color()`

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
        fl.basename = nodename  # this is a good example of how we can
        #                         use basename to make it easy to deal
        #                         with multiple datasets
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
            s = prscr.clock_correct(s)
        # TODO ☐: why do you return dropped?
        s, dropped_dims = s.squeeze(return_dropped=True)
        s, ax_last = prscr.rough_table_of_integrals(s, fl=fl)
        # Included signal averaging in rough_table_of_integrals
        Mi, R1, vd = sympy.symbols("M_inf R_1 vd", real=True)
        psd.logger.debug(psd.strm("acq keys", s.get_prop("acq_params")))
        W = (
            s.get_prop("acq_params")["FIR_rep"] * 1e-6
            + s.get_prop("acq_params")["acq_time_ms"] * 1e-3
        )
        s = psd.lmfitdata(s)
        s.functional_form = Mi * (
            1 - (2 - sympy.exp(-W * R1)) * sympy.exp(-vd * R1)
        )
        s.set_guess(
            M_inf=dict(
                value=s.max().item(),
                min=0.1 * s.max().item(),
                max=1.5 * s.max().item(),
            ),
            R_1=dict(value=0.8, min=0.01, max=100),
        )
        s.fit()
        print(s.get_units("vd"))
        s_fit = s.eval(200)
        print(s_fit.get_units("vd"))
        psd.plot(s_fit, ax=ax_last, alpha=0.5)  # here, we plot the fit
        #                                         together with the
        #                                         table of integrals.
        if plot_fit:  # JF has not reviewed this -- needs to be re-written
            #       consistently w/ above.  Stuff that's not used can
            #       just be removed
            T1 = 1.0 / s.output("R_1")
            Mi = s.output("M_inf")
            T1_list.append(T1)
            s_fit.set_units("vd", "s")
            fit = s.eval(100)
            fit.set_plot_color(s_fit.get_plot_color())
            fl.basename = None  # because we want the following plot to
            #                    show up together
            fl.next("IR fit - normalized")
            fl.plot(s / Mi, "o", label=nodename)
            fl.plot(
                fit / Mi,
                ls="-",
                alpha=0.5,
                label="fit for %s" % nodename,
            )
            ax = plt.gca()
# I'm not printing anything for 'T1 = ?' as desired in the list of goals, what
# should I be printing? T1 at s.max()?
# TODO ☐: pyspecdata gives an example of how to show the fit equation on
#         the plot (fit_fake_data.py).  I would do that on the bottom
#         right plot of the table of integrals.
