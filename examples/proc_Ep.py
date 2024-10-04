"""
Process Enhancement experiment 
====================================================
Opens .h5 results file, uses rough_table_of_integrals() to roughly process
dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
import os
import datetime
import matplotlib.pyplot as plt

data_target = os.path.normpath(psd.getDATADIR("WK_processed_data"))

with psd.figlist_var() as fl:
    thisfile, exptype, nodename, log_nodename = (
        "240924_13p5mM_TEMPOL_ODNP_1.h5",
        "ODNP_NMR_comp/ODNP",
        "ODNP",
        "log",
    )
    # TODO ☐: it complains about the following, as it should,
    #         because it's not an nddata.  This is why
    #         I was saying you will want to pass the
    #         filename and node name to your function.
    s_log = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=log_nodename,
    )
    s_log, poweraxis = prscr.convert_to_power(s_log, fl=fl)
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        lookup=prscr.lookup_table,
    )
    orig_axis = s["indirect"] # let's save this so we
    #                           can pass it to the log
    s["indirect"] = s["indirect"]["start_times"]
    s.set_units("indirect", "s")
    s, _ = prscr.rough_table_of_integrals(s, fl=fl)
    s.set_error(s["indirect", 0].item() * 0.01)  # We are not calculating the
    #                                              errors in rough table of
    #                                              integrals, so just make up a
    #                                              reasonable sized random
    #                                              number so that I can see the
    #                                              relative errors!
    s /= s["indirect", 0:1]
    fl.next("normalized $E(p(t))$")
    s["indirect"] -= s["indirect"][0]
    fl.plot(s, "o")
    # {{{ this is just matplotlib time formatting
    ax = plt.gca()
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(
            lambda x, _: str(datetime.timedelta(seconds=x))
        )
    )
    # }}}
