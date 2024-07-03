import numpy as np
import pyspecdata as psd
from pyspecProcScripts import lookup_table, select_pathway

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

with psd.figlist_var() as fl:
    for thisfile, exp_type, nodename, this_basename, manystep_cpmg, thislabel in [
        (
            "240702_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_13",
            True,
            False,
            "CPMG simple cyc SW = 16",
        ),
        (
            "240702_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_19",
            True,
            True,
            "CPMG large cyc SW = 10.0",
        ),
        (
            "240702_13p5mM_TEMPOL_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_12",
            True,
            False,
            "echo large cyc SW = 10.0",
        ),
    ]:
        fl.this_basename = no_power

        thisd = psd.find_file(
            thisfile, exp_type=exp_type, expno=nodename, lookup=lookup_table
        )
        forplot = thisd.C.sum("nScans").ift("t2").run(abs)
        thisd.squeeze()
        fl.next("raw data for %s" % thislabel)
        fl.image(thisd, interpolation="none")
        fl.next("abs of raw data for %s, signal average" % thislabel)
        if "nEcho" in forplot.dimlabels:
            forplot.smoosh(["nEcho", "t2"], r"nEcho $\otimes$ t2")
        fl.image(forplot, interpolation="none")
        thisd.ift("t2")
        fl.next("abs(t domain) comparison")
        if manystep_cpmg:
            thisd = echo_interleave(thisd, "ph1")
            thisd.set_prop("coherence_pathway", {"ph2": -2, "ph_overall": -1})
        thisd = select_pathway(thisd, thisd.get_prop("coherence_pathway"))
        if "nEcho" in thisd.dimlabels:
            thisd.smoosh(["nEcho", "t2"], "t2")
            acq = thisd.get_prop("acq_params")
            echo_time = 1e-6 * 2 * (acq["tau_us"] + acq["p90_us"])
            thisd["t2"] = (thisd["t2"]["nEcho"]) * echo_time + thisd["t2"]["t2"]
        if 'nScans' in thisd.dimlabels:
            thisd.mean('nScans')
        fl.plot(abs(thisd), 'o',label=thislabel)
