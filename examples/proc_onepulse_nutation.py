import pyspecdata as psd
import pyspecProcScripts as prscr
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import sympy as sp
from sympy import symbols

signal_range = (-250, 0)


def clock_correct(s, axis_along, direct="t2", max_cyc=0.5):
    for_correct = s.C
    Delta = np.diff(s[axis_along][np.r_[0, -1]]).item()
    correction_axis = psd.nddata(np.r_[-0.5:0.5:300j] * max_cyc / Delta, "correction")
    for_correct = for_correct * np.exp(
        -1j * 2 * pi * correction_axis * for_correct.fromaxis(axis_along)
    )
    # {{{ determine the best sign flip for each correction
    for j in range(for_correct.shape["correction"]):
        thesign = prscr.determine_sign(for_correct["correction", j])
        for_correct["correction", j] *= thesign
    # }}}
    for_correct.sum(direct)
    return for_correct.sum(axis_along).run(abs).argmax("correction").item()


with psd.figlist_var() as fl:
    for searchstr, postproc_type, nodename in [
        ["240311_27mM_TEMPOL_nutation_1.h5", "spincore_nutation_v2", "nutation"]
    ]:
        s = psd.find_file(
            searchstr,
            exp_type="ODNP_NMR_comp/nutation",
            expno=nodename,
            lookup=prscr.lookup_table,
            postproc=postproc_type,
        )
        s["indirect"] *= 1e-6  # convert to s
        s.set_prop("coherence_pathway", {"ph1": 1})
        s.rename("indirect", "p_90")
        s.mean("nScans")
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        fl.next("Raw Data")
        s.reorder(["ph1", "p_90", "t2"])
        fig.suptitle("Nutation Processing")
        fl.image(s["t2":(-500, 500)].C, ax=ax1)
        ax1.set_title("Raw Data")
        mysgn = prscr.determine_sign(prscr.select_pathway(s["t2":signal_range]))
        corr = clock_correct(prscr.select_pathway(s), "p_90")
        s /= prscr.zeroth_order_ph(
            prscr.select_pathway(s["t2":signal_range].C.sum("t2"))
        )
        my_sign = prscr.determine_sign(prscr.select_pathway(s["t2":signal_range]))
        s *= my_sign
        s = prscr.fig_from_echo(s, s.get_prop("coherence_pathway"))
        for j in range(len(s.getaxis("p_90"))):
            s["p_90", j] /= prscr.zeroth_order_ph(
                prscr.select_pathway(s["p_90", j]["t2":signal_range]).C.sum("t2")
            )
        s.ift("t2")
        s = s["t2":(0, None)]
        s *= 2
        s["t2":0] *= 0.5
        s.ft("t2")
        fl.image(s["t2":disp_range].C, ax=ax2)
        ax2.set_title("Phased and FID sliced")
        s *= my_sign
        fl.image(s, ax=ax3)
        ax3.set_title("Sign Flipped Back")
