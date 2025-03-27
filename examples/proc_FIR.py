"""
Process FIR experiment 
======================
Opens .h5 results file, uses :func:`rough_table_of_integrals` to roughly
process
dataset including generating a table of integrals
"""

import pyspecProcScripts as prscr
import pyspecdata as psd
from pyspecdata import Q_
import sympy
import matplotlib.pyplot as plt
import numpy as np
import re
import h5py

plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2


clock_correction = True
plot_fit = True
thisfile, exptype, post_proc, lookup = (
    "240924_13p5mM_TEMPOL_ODNP_1.h5",
    "ODNP_NMR_comp/ODNP",
    "spincore_IR_v3",
    prscr.lookup_table,
)
filename = psd.search_filename(
    re.escape(thisfile),
    exp_type=exptype,
    unique=True,
)
with h5py.File(filename, mode="r") as fp:
    R1nodenames = [j for j in fp.keys() if "FIR" in j]
# Because we are going ot want to get both R1 fit values as well as the
# associated errors, we collect the results in an nddata rather than
# just e.g. a list
R1data = psd.ndshape([("power", len(R1nodenames))]).alloc(dtype=np.float64)

with psd.figlist_var() as fl:
    for j, nodename in enumerate(R1nodenames):
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
        s.set_plot_color_next()
        indirect = "vd"
        direct = "t2"
        if clock_correction:
            s = prscr.clock_correct(s)
        s = s.squeeze()
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
        prefactor_scaling = s.div_units("vd", "s")
        length_of_decay = s["vd"].max()
        int_at_end_of_decay = s["vd":length_of_decay].item()
        s.set_guess(
            # here, rather than trying a max or min, we pull the last
            # point along the indirect (since it can vary in sign)
            M_inf=dict(
                value=int_at_end_of_decay,
                min=0.1 * int_at_end_of_decay,
                max=5 * int_at_end_of_decay,
            ),
            R_1=dict(
                value=5 / length_of_decay,
                min=0.01 / length_of_decay,
                max=100 / length_of_decay,
            ),
        )
        s.fit()
        s_fit = s.eval(200)
        psd.plot(s_fit, ax=ax_last, alpha=0.5)  # here, we plot the fit
        #                                         together with the
        #                                         table of integrals.
        R1 = s.output("R_1")
        ax_last.text(
            0.5,
            0.5,
            (
                "%s RESULT: %s\n$R_1=%#0.3g\\;" % (nodename, s.latex(), R1)
                + f"{1/Q_(s.get_units('vd')):~L}$"
            ),
            ha="center",
            va="center",
            color=s_fit.get_plot_color(),
            transform=ax_last.transAxes,
        )
        if plot_fit:  # JF has not reviewed this -- needs to be re-written
            #       consistently w/ above.  Stuff that's not used can
            #       just be removed
            R1data["power", j] = R1
            Mi = s.output("M_inf")
            fl.basename = None  # because we want the following plot to
            #                    show up together
            fl.next("IR fit - normalized")
            fl.plot(s / Mi, "o", label=nodename)
            fl.plot(
                s_fit / Mi,
                ls="-",
                alpha=0.5,
                label="fit for %s" % nodename,
            )
            ax = plt.gca()
# I'm not printing anything for 'T1 = ?' as desired in the list of goals, what
# should I be printing? T1 at s.max()?
