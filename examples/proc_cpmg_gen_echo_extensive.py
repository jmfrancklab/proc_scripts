import numpy as np
import pyspecdata as psd
import matplotlib.pylab as plt
from itertools import cycle
from pyspecProcScripts import lookup_table, select_pathway


def size_of_coherences(thisd, fl):
    forplot = thisd.C
    if "nEcho" in forplot.dimlabels:
        forplot.smoosh(["nEcho", "t2"], r"nEcho $\otimes$ t2")
        if "nScans" in forplot.dimlabels:
            lastph = [j for j in forplot.dimlabels if j.startswith("ph")][-1]
            forplot.smoosh([lastph, "nScans"], r"last $\Delta p$ $\otimes$ nScans")
        forplot.reorder([r"nEcho $\otimes$ t2"], first=False)
    fl.image(forplot)


def echo_interleave(d, phcycdim):
    "interleave even and odd echoes coming from phcycdim"
    assert d.get_ft_prop(phcycdim)
    retval = d.shape
    retval.pop(phcycdim)
    retval = retval.alloc(dtype=np.complex128, format=None)
    retval["nEcho", 0::2] = d["ph1", +1]["nEcho", 0::2]
    retval["nEcho", 1::2] = d["ph1", -1]["nEcho", 1::2]
    retval.copy_axes(d).copy_props(d)
    return retval


def smoosh_direct_domain(thisd):
    thisd.smoosh(["nEcho", "t2"], "t2")
    acq = thisd.get_prop("acq_params")
    echo_time = 1e-6 * 2 * (acq["tau_us"] + acq["p90_us"])
    thisd["t2"] = (thisd["t2"]["nEcho"]) * echo_time + thisd["t2"]["t2"]
    return thisd


colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
plt.rcParams.update(
    {
        "errorbar.capsize": 2,
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
    }
)
colorcyc = cycle(colorcyc_list)


with psd.figlist_var() as fl:
    for (
        thisfile,
        exp_type,
        nodename,
        manystep_cpmg,
        thislabel,
        thisbasename,
    ) in [
        (
            "240702_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_19",
            True,
            "cpmg19 no power large cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_1",
            True,
            "pm cpmg1 no power large cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_2",
            True,
            "pm cpmg2 no power large cyc SW = 10.0",
            "no power",
        ),
    ]:
        fl.basename = thisbasename
        thisd = psd.find_file(
            thisfile, exp_type=exp_type, expno=nodename, lookup=lookup_table
        )
        thisd.squeeze()
        fl.next("raw data for %s" % thislabel)
        fl.image(thisd, interpolation="auto")
        fl.next("abs of raw data for %s, signal average" % thislabel)
        thisd.ift("t2")
        size_of_coherences(thisd, fl)
        if manystep_cpmg:
            proposed_cycle = thisd.C.ift(["ph1", "ph2", "ph_overall"])
            proposed_cycle = proposed_cycle["ph_overall", 0::2]["ph2", 0::2]
            proposed_cycle.ft(["ph1", "ph2", "ph_overall"])
            fl.next("signal from proposed cycle")
            size_of_coherences(proposed_cycle, fl)
            fl.next("reduced phcyc with interleave")
            temp = echo_interleave(proposed_cycle, "ph1")
            temp = smoosh_direct_domain(temp)
            fl.plot(abs(temp["ph2", -2]["ph_overall", -1]), "o")
            fl.next("signal from CPMG-only")
            proposed_cycle = proposed_cycle["ph1", 1] - proposed_cycle["ph1", -1]
            proposed_cycle = smoosh_direct_domain(proposed_cycle)
            fl.plot(abs(proposed_cycle["ph2", -2]["ph_overall", -1]), "o")
            thisd = echo_interleave(thisd, "ph1")
            thisd.set_prop("coherence_pathway", {"ph2": -2, "ph_overall": -1})
        thisd = select_pathway(thisd, thisd.get_prop("coherence_pathway"))
        if "nEcho" in thisd.dimlabels:
            thisd = smoosh_direct_domain(thisd)
        fl.next("plot each scan separately")
        # c = next(colorcyc)
        kwargs = {}
        for j in range(thisd.shape["nScans"]):
            fl.plot(abs(thisd["nScans", j]), "o", label=f"{thislabel} scan {j}")
        fl.next("abs(t domain) comparison")
        if "nScans" in thisd.dimlabels:
            thisd.mean("nScans")
        fl.plot(abs(thisd), "o", label=thislabel)
