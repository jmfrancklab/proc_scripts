"""
Process FIR experiment 
====================================================
Opens .h5 results file, uses rough_table_of_integrals() to roughly process
dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
import numpy as np
import sympy
import matplotlib.pyplot as plt
from itertools import cycle

color_cycle = cycle(
    [
        "red",
        "orange",
        "yellow",
        "green",
        "cyan",
        "blue",
        "purple",
        "magenta",
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
)
T1_list = []

with psd.figlist_var() as fl:
    thisfile, exptype, post_proc, lookup = (
        "240924_13p5mM_TEMPOL_ODNP_1.h5",
        "ODNP_NMR_comp/ODNP",
        "spincore_IR_v3",
        prscr.lookup_table,
    )
    for nodename in [
        # "FIR_noPower",
        "FIR_34dBm",
    ]:
        print(nodename)
        s = psd.find_file(
            thisfile,
            exp_type=exptype,
            expno=nodename,
            postproc=post_proc,
            lookup=prscr.lookup_table,
        )
        s.set_units("vd", "s")
        clock_correction = True
        indirect = "vd"
        direct = "t2"
        if clock_correction:
            # {{{ clock correction
            clock_corr = psd.nddata(np.linspace(-3, 3, 2500), "clock_corr")
            assert s.get_ft_prop(direct)
            if fl is not None:
                fl.next("before clock correction")
                fl.image(s)
            s_clock = (
                prscr.select_pathway(s, s.get_prop("coherence_pathway"))
                .mean("nScans")
                .sum(direct)
            )
            phase_dims = [j for j in s.dimlabels if j.startswith("ph")]
            s.ift(phase_dims)
            min_index = abs(s_clock).argmin(indirect, raw_index=True).item()
            s_clock *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
            s_clock[indirect, : min_index + 1] *= -1
            s_clock.sum(indirect).run(abs)
            if fl is not None:
                fl.next("clock correction")
                fl.plot(s_clock, ".", alpha=0.7)
            clock_corr = s_clock.argmax("clock_corr").item()
            plt.axvline(x=clock_corr, alpha=0.5, color="r")
            s *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
            s.ft(phase_dims)
            if fl is not None:
                fl.next("after auto-clock correction")
                fl.image(s)
            # }}}
        np.squeeze(s)
        s, _ = prscr.rough_table_of_integrals(s, fl=fl)
        if "nScans" in s.dimlabels:
            s_int = s.C.mean("nScans")
        else:
            s_int = s.C
        print(s_int)
        thiscolor = next(color_cycle)
        Mi, R1, vd = sympy.symbols("M_inf R_1 vd", real=True)
        psd.logger.debug(psd.strm("acq keys", s.get_prop("acq_params")))
        W = (
            s.get_prop("acq_params")["FIR_rep"] * 1e-6
            + s.get_prop("acq_params")["acq_time_ms"] * 1e-3
        )
        functional_form = Mi * (
            1 - (2 - sympy.exp(-W * R1)) * sympy.exp(-vd * R1)
        )
        IR_data = psd.nddata(s_int.C.data, ["vd"])
        IR_data.setaxis("vd", s_int.getaxis("vd"))
        f = psd.lmfitdata(IR_data)
        f.functional_form = functional_form
        f.set_guess(
            M_inf=dict(value=3.9e4, min=0, max=2e7),
            R_1=dict(value=0.8, min=0.01, max=100),
        )
        fl.next("IR fit - before norm")
        fl.plot(s_int, "o", color=thiscolor, label="%s" % nodename)
        f.fit()
        T1 = 1.0 / f.output("R_1")
        Mi = f.output("M_inf")
        T1_list.append(T1)
        f.set_units("vd", "s")
        fit = f.eval(100)
        fl.plot(
            fit,
            ls="-",
            color=thiscolor,
            alpha=0.5,
            label="fit for %s" % nodename,
        )
        fl.next("IR fit - before norm - %s" % nodename)
        fl.plot(s_int, "o", color=thiscolor, label="%s" % nodename)
        fl.plot(
            fit,
            ls="-",
            color=thiscolor,
            alpha=0.5,
            label="fit for %s" % nodename,
        )
        ax = plt.gca()
