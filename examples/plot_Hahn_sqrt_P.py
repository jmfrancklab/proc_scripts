import pyspecdata as psd
import pyspecProcScripts as prscr
from itertools import cycle
import matplotlib.pyplot as plt
import numpy as np

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)
dg_per_V_3p9 = 1.257e11
fig, ax_list = plt.subplots(2, 2, sharey=True)
with psd.figlist_var() as fl:
    for searchstr, label, nodename, expected, ax_set in [
        [
            "240805_amp1_27mM_TEMPOL_echo.h5",
            "amplitude = 1",
            "echo_4",
            6.61,
            ax_list[0][0],
        ],
        [
            "240805_amp0p2_27mM_TEMPOL_echo.h5",
            "amplitude = 0.2",
            "echo_4",
            7.9,
            ax_list[0][1],
        ],
        [
            "240805_amp0p1_27mM_TEMPOL_echo.h5",  # optimal beta
            "amplitude = 0.1",
            "echo_4",
            6,
            ax_list[1][0],
        ],
        [
            "240805_amp0p05_27mM_TEMPOL_echo.h5",
            "amplitude = 0.05",
            "echo_3",
            6,
            ax_list[1][1],
        ],
    ]:
        s = psd.find_file(
            searchstr,
            exp_type="ODNP_NMR_comp/Echoes",
            expno=nodename,
            lookup=prscr.lookup_table,
        )
        if "nScans" in s.dimlabels:
            s.mean("nScans")
        s /= dg_per_V_3p9
        s.ift("t2")
        s /= np.sqrt(50)
        s.ft("t2")
        # {{{ phasing
        ph0 = s["t2":(-500, 500)].C.sum("t2")
        ph0 /= abs(ph0)
        s /= ph0
        # }}}
        s.ift("t2")
        # {{{ convert to sqrt(fW)
        s /= np.sqrt(1e-15)
        # }}}
        s = abs(s)
        fl.next("Echo amplitude vs artefacts")  # %s'%label,legend=True)
        fl.plot(
            prscr.select_pathway(s["t2":(0, 150e-3)], s.get_prop("coherence_pathway")),
            alpha=0.5,
            label="absolute signal",
            ax=ax_set,
        )
        ax_set.set_title("%s" % label)
        s["ph1", 1] = 0
        s.sum("ph1")
        s.set_units("t2", "s")
        fl.plot(
            s["t2":(0, 150e-3)],
            alpha=0.5,
            label="off coherence pathways summed",
            ax=ax_set,
            human_units=False,
        )
        plt.ylabel("$\\sqrt{P}$ / $\\mathrm{\\sqrt{fW}}$")
        plt.xlabel("$\\tau$ / ms")
        plt.legend()
        ax_set.grid(which="both")
        plt.legend()
        ax_set.axhline(expected)
