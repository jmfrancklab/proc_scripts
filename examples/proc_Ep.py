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
import re

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
    # TODO ☐: make into a loop that includes this -- I only didn't want
    #         to change the indent b/c I wanted you to be able to see the diff
    # thisfile, exptype, nodename = (
    #    "240924_13p5mM_TEMPOL_ODNP_1.h5",
    #    "ODNP_NMR_comp/ODNP",
    #    "ODNP",
    # )
    thisfile, thisexptype, nodename = (
        "260429_hydroxytempo_ODNP_2.h5",
        "B27/ODNP",
        "ODNP",
    )
    s = psd.find_file(
        thisfile,
        exp_type=thisexptype,
        expno=nodename,
        lookup=prscr.lookup_table,
    )
    orig_axis = s["indirect"]  # let's save this so we
    #                           can pass it to the log
    orig_axis_error = s.get_error("indirect")
    s["indirect"] = s["indirect"][
        "time"
    ]  # we need to do this so that the rough table of
    s.set_error("indirect", orig_axis_error["time"])
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
    s.set_error("indirect", orig_axis_error)
    m = re.search(r".*ODNP.*v([0-9]+)$", s.get_prop("postproc_type"))
    if m:
        vernum = int(m.group(1))
    else:
        raise IOError("What the heck type of postproc type is that!!")
    if s.get_prop("log") is None:
        # the log has not changed, just its location in the HDF file.
        # so, once we attach it, the following should work!
        s = prscr.attach_log_data_from_file(s, thisfile, thisexptype)
    if vernum < 6:
        # the next line is done already if we are 6 or higher
        s = prscr.generate_coordinates_from_log(s, fl=fl)
    # {{{ our standard philosophy for postproc is that we do not destroy
    #    info -- so the following, which destroys info about everything
    #    but power, belongs here.
    s.set_error("indirect", s.get_error("indirect")["power"])
    s["indirect"] = s["indirect"]["power"]
    s.set_units("indirect", "W").rename("indirect", "power")
    # }}}
    fl.next("normalized $E(p)$")
    fl.plot(s, "o")
