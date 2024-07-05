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
    for thisfile, exp_type, nodename, manystep_cpmg, thislabel, thisbasename in [
        # {{{ Tested Controls
        (
            "240704_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_1",
            True,
            "CPMG no power large cyc SW = 10.0",
            "control",
        ),
        (
            "240704_13p5mM_TEMPOL_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_1",
            False,
            "echo no power large cyc SW = 10.0",
            "control",
        ),
        (
            "240704_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_2",
            False,
            "CPMG no power simple cyc SW = 10.0",
            "control",
        ),
        # }}}
        # {{{ No power
        (
            "240702_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_19",
            True,
            "CPMG no power large cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_12",
            False,
            "echo no power large cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_1",
            True,
            "CPMG no power simple cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_2",
            True,
            "CPMG no power large cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_2",
            False,
            "echo no power large cyc SW = 10.0",
            "no power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_3",
            False,
            "echo no power large cyc SW = 10.0",
            "no power",
        ),
        # }}}
        # {{{ 30 dBm
        (
            "240702_13p5mM_TEMPOL_pm_30dB_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_1",
            True,
            "CPMG 30dBm large cyc SW = 10.0",
            "30 dBm uw power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_30dB_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_1",
            False,
            "echo 30dBm large cyc SW = 10.0",
            "30 dBm uw power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_30dB_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_2",
            True,
            "CPMG 30dBm large cyc SW = 10.0",
            "30 dBm uw power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_30dB_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_2",
            False,
            "echo 30dBm large cyc SW = 10.0",
            "30 dBm uw power",
        ),
        # }}}
        # {{{ 34 dBm
        (
            "240702_13p5mM_TEMPOL_pm_34dB_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_1",
            True,
            "CPMG 34dBm large cyc SW = 10.0",
            "34 dBm uw power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_34dB_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_1",
            False,
            "echo 34 dBm large cyc SW = 10.0",
            "34 dBm uw power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_34dB_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_2",
            True,
            "CPMG 34 dBm large cyc SW = 10.0",
            "34 dBm uw power",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_34dB_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_2",
            False,
            "echo 34 dBm large cyc SW = 10.0",
            "34 dBm uw power",
        ),
        # }}}
        #{{{ New Probe vs Old
        (
            "240705_13p5mM_TEMPOL_balancedprobe_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_11",
            False,
            "echo balanced probe 40W",
            "compare_probes",
        ),
        (
            "240705_13p5mM_TEMPOL_VTUprobe_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_14",
            False,
            "echo VTU probe 40W",
            "compare_probes",
        ),
        #}}}
    ]:
        fl.basename = thisbasename
        thisd = psd.find_file(
            thisfile, exp_type=exp_type, expno=nodename, lookup=lookup_table
        )
        thisd.squeeze()
        fl.next("raw data for %s" % thislabel)
        fl.image(thisd, interpolation="auto")
        fl.next("abs of raw data for %s, signal average" % thislabel)
        forplot = thisd.C.sum("nScans").ift("t2").run(abs)
        if "nEcho" in forplot.dimlabels:
            forplot.smoosh(["nEcho", "t2"], r"nEcho $\otimes$ t2")
        fl.image(forplot, interpolation="auto")
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
        if "nScans" in thisd.dimlabels:
            thisd.mean("nScans")
        fl.plot(abs(thisd), "o", label=thislabel)
