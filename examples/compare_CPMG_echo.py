import numpy as np
from pyspecdata import find_file, figlist_var
from pyspecProcScripts import lookup_table, select_pathway

with figlist_var() as fl:
    for thisfile, exp_type, nodename, complex_cpmg, thislabel in [
        (
            "240702_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_13",
            False,
            "CPMG simple cyc",
        ),
        (
            "240702_13p5mM_TEMPOL_CPMG.h5",
            "ODNP_NMR_comp/CPMG",
            "CPMG_16",
            True,
            "CPMG large cyc",
        ),
        (
            "240702_13p5mM_TEMPOL_echo.h5",
            "ODNP_NMR_comp/Echoes",
            "echo_7",
            False,
            "echo large cyc",
        ),
    ]:
        thisd = find_file(
            thisfile, exp_type=exp_type, expno=nodename, lookup=lookup_table
        )
        thisd.squeeze()
        fl.next("raw data for %s" % thislabel)
        rows = np.prod([thisd.shape[j] for j in thisd.dimlabels[:-1]])
        if rows < 500:
            fl.image(thisd)
        else:
            fl.image(thisd, interpolation="bilinear")
        thisd.ift("t2")
        fl.next("abs(t domain) comparison")
        if complex_cpmg:
            odd = select_pathway(thisd,thisd.get_prop('coherence_pathway'))
            even_pw = {'ph1':-1,'ph2':-2,'ph_overall':-1}
            even = select_pathway(thisd.C,even_pw)
            thisd = odd + even
        else:
            thisd = select_pathway(thisd, thisd.get_prop("coherence_pathway"))
        if "nEcho" in thisd.dimlabels:
            thisd.smoosh(["nEcho", "t2"], "t2")
            acq = thisd.get_prop("acq_params")
            echo_time = 1e-6 * 2 * (acq["tau_us"] + acq["p90_us"])
            thisd["t2"] = (thisd["t2"]["nEcho"]) * echo_time + thisd["t2"]["t2"]
        fl.plot(abs(thisd), label=thislabel)
