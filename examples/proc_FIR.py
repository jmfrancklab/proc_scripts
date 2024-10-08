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
# I'm using the color cycle for more plots than just those made by
# rough_table_of_integrals, I don't have this iterating correctly though
prop_cycle = plt.rcParams["axes.prop_cycle"]
color_cycle = prop_cycle.by_key()["color"]
thiscolor = next(iter(color_cycle))

# should I be using something like this instead?
# colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# colorcyc = cycle(colorcyc_list)


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
            M_inf=dict(value=s.max().item(), min=0, max=1.5 * s.max().item()),
            R_1=dict(value=0.8, min=0.01, max=100),
        )
        s.fit()
        s_fit = s.eval(200)
        psd.plot(s_fit, ax=ax_last, alpha=0.5)  # here, we plot the fit
        #                                       together with the table of
        #                                       integrals.
        if plot_fit:  # JF has not reviewed this -- needs to be re-written
            #       consistently w/ above.  Stuff that's not used can
            #       just be removed
            T1 = 1.0 / s.output("R_1")
            Mi = s.output("M_inf")
            T1_list.append(T1)
            s_fit.set_units("vd", "s")
            fit = s.eval(100)
            fit.set_plot_color(s_fit.get_plot_color())
            fl.next("IR fit - normalized - %s" % nodename)
            fl.plot(s / Mi, "o", color=thiscolor, label=nodename)
            fl.plot(
                fit / Mi,
                ls="-",
                color=thiscolor,
                alpha=0.5,
                label="fit for %s" % nodename,
            )
            ax = plt.gca()
# I'm not printing anything for 'T1 = ?' as desired in the list of goals, what
# should I be printing? T1 at s.max()?
