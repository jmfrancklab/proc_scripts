import numpy as np
from pyspecdata import find_file, figlist_var
from pyspecProcScripts import lookup_table, select_pathway

with figlist_var() as fl:
    for thisfile, exp_type, nodename, complex_cpmg, thislabel in [
        #(
        #    "240702_13p5mM_TEMPOL_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_13",
        #    False,
        #    "CPMG simple cyc",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_17",
        #    True,
        #    "CPMG large cyc SW = 2.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_19",
        #    True,
        #    "CPMG large cyc SW = 10",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_10",
        #    False,
        #    "echo large cyc SW = 2.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_12",
        #    False,
        #    "echo large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_1",
        #    False,
        #    "pm echo large cyc SW = 3.9",
        #),
        #{{{ B12 on, no power
        #(
        #    "240702_13p5mM_TEMPOL_pm_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_2",
        #    False,
        #    "pm echo large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_3",
        #    False,
        #    "pm echo large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_1",
        #    True,
        #    "pm cpmg large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_2",
        #    True,
        #    "pm cpmg large cyc SW = 10.0",
        #),
        #}}}
        #{{{ 34 dB uw
        #(
        #    "240702_13p5mM_TEMPOL_pm_34dB_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_1",
        #    False,
        #    "pm echo 34 dB large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_34dB_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_1",
        #    True,
        #    "pm cpmg 34 dB large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_34dB_echo.h5",
        #    "ODNP_NMR_comp/Echoes",
        #    "echo_2",
        #    False,
        #    "pm echo 34 dB large cyc SW = 10.0",
        #),
        #(
        #    "240702_13p5mM_TEMPOL_pm_34dB_CPMG.h5",
        #    "ODNP_NMR_comp/CPMG",
        #    "CPMG_2",
        #    True,
        #    "pm cpmg 34 dB large cyc SW = 10.0",
        #),
        #}}}
        #{{{ 30 dB uw
        (
            "240702_13p5mM_TEMPOL_pm_30dB_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_1",
            False,
            "pm echo1 30 dB large cyc SW = 10.0",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_30dB_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_1",
            True,
            "pm cpmg1 30 dB large cyc SW = 10.0",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_30dB_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_2",
            False,
            "pm echo2 30 dB large cyc SW = 10.0",
        ),
        (
            "240702_13p5mM_TEMPOL_pm_30dB_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_2",
            True,
            "pm cpmg2 30 dB large cyc SW = 10.0",
        ),
  ]:
        thisd = find_file(
            thisfile, exp_type=exp_type, expno=nodename, lookup=lookup_table
        )
        thisd.squeeze()
        fl.next("raw data for %s" % thislabel)
        fl.image(thisd, interpolation='auto')
        fl.next("raw data for %s, slice overall and signal average" % thislabel)
        forplot = thisd['ph_overall',-1].sum('nScans').ift('t2')
        if 'nEcho' in forplot.dimlabels:
            forplot.reorder('nEcho', first=False)
        fl.image(forplot)
        thisd.ift("t2")
        fl.next("abs(t domain) comparison")
        if complex_cpmg:
            odd_pw = {'ph1':1,'ph2':-2,'ph_overall':-1}
            odd = select_pathway(thisd,odd_pw)
            even_pw = {'ph1':-1,'ph2':+2,'ph_overall':-1}
            even = select_pathway(thisd.C,even_pw)
            # {{{ interleave the echoes
            thisd = odd
            # "even" starts w/ 1 b/c that's the second echo
            thisd['nEcho',1::2] = even['nEcho',1::2]
            # }}}
        else:
            thisd = select_pathway(thisd, thisd.get_prop("coherence_pathway"))
        if "nEcho" in thisd.dimlabels:
            thisd.smoosh(["nEcho", "t2"], "t2")
            acq = thisd.get_prop("acq_params")
            echo_time = 1e-6 * 2 * (acq["tau_us"] + acq["p90_us"])
            thisd["t2"] = (thisd["t2"]["nEcho"]) * echo_time + thisd["t2"]["t2"]
        if 'nScans' in thisd.dimlabels:
            thisd.mean('nScans')
        fl.plot(abs(thisd), 'o',label=thislabel)
