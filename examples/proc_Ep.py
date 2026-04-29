"""
Process Enhancement experiment
====================================================
Opens .h5 results file, uses rough_table_of_integrals() to roughly process
dataset including generating a table of integrals
"""

import pyspecProcScripts as prscr
import pyspecdata as psd
import datetime
import matplotlib.pyplot as plt

plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2
plt.rcParams.update(
    {
        "errorbar.capsize": 2,
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.figsize": (6, 5),
    }
)

with psd.figlist_var() as fl:
    thisfile, exptype, nodename = (
        "260429_hydroxytempo_ODNP_1.h5",
        "B27/ODNP",
        "ODNP",
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        lookup=prscr.lookup_table,
    )
    orig_axis = s["indirect"]  # let's save this so we
    #                           can pass it to the log
    s["indirect"] = (
        s["indirect"]["start_times"] - s["indirect"]["start_times"][0]
    )
    s.set_units("indirect", "s")
    s, _ = prscr.rough_table_of_integrals(s, fl=fl)
    assert psd.det_unit_prefactor(s.get_units("indirect")) == 0
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
        plt.FuncFormatter(lambda x, _: str(datetime.timedelta(seconds=x)))
    )
    # }}}
    s["indirect"] = orig_axis
    post_proc_type = s.get_prop("postproc_type")
    if post_proc_type == "spincore_ODNP_v6":
        s = prscr.generate_coordinates_from_log(s, thisfile, exptype, fl=fl)
        indirect_error = s.get_error("indirect")
        s.setaxis("indirect", s["indirect"]["power"])
        s.set_error("indirect", indirect_error["power"])
        s.set_units("indirect", "dBm").rename("indirect", "power")
    else:
        s = prscr.generate_power_coordinates_from_log(
            s,
            thisfile,
            exptype,
            fl=fl,
        )
    fl.next("normalized $E(p)$")
    fl.plot(s, "o")
