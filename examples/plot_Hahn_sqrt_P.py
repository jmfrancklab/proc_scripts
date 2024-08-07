import pyspecdata as psd
import pyspecProcScripts as prscr
from itertools import cycle
import matplotlib.pyplot as plt
import numpy as np

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)
dg_per_V_3p9 = 1.257e11
with psd.figlist_var() as fl:
    for searchstr, label, nodename in [
        ["240805_amp0p05_27mM_TEMPOL_echo.h5", "amplitude = 0.05", "echo_3"],
    ]:
        s = psd.find_file(
            searchstr,
            exp_type="ODNP_NMR_comp/Echoes",
            expno=nodename,
            lookup=prscr.lookup_table,
        )
        if "nScans" in s.dimlabels:
            s.mean("nScans")
        s.set_units("t2", "s")
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
        fl.next("Echo admplitude 0.05")  # %s'%label,legend=True)
        fl.plot(
            prscr.select_pathway(s["t2":(0, 150e-3)], s.get_prop("coherence_pathway")),
            alpha=0.5,
            label="absolute signal %s" % label,
        )  # for j in range(s.shape['ph1']):
        plt.legend()
    plt.xlabel("$\\tau$ / ms")
    plt.ylabel("$\\sqrt{P}$ / $\\mathrm{\\sqrt{fW}}$")
    plt.legend()
    psd.gridandtick(plt.gca())
