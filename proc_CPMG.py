from pyspecdata import *
from scipy.optimize import basinhopping
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from pyspecProcScripts.fitting import decay
import sympy as sp

fl = fl_mod()
logger = init_logging("info")
for (
    searchstr,
    exp_type,
    nodename,
    postproc,
    label_str,
    f_range,
) in [
    (
        "w8_200917",
        "test_equip",
        6,
        "ag_CPMG_strob",
        "water loading 8",
        (-500, 500),
    ),
    # ('200221_CPMG_TEMPOLgel_2p9_1','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('freeSL_201007','test_equip',3,'ag_CPMG_strob','free SL',(-200,200)),
    # ('200305_CPMG_3p5_2','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('200305_CPMG_3p6_2','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('200305_CPMG_3p7_2','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('200305_CPMG_3p7_3','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('200305_CPMG_3p8_2','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('200305_CPMG_3p9_2','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
    # ('200305_CPMG_4p0_1','test_equip','signal','spincore_CPMG_v1',
    #    'deadtime=5',(-500,500)),
]:
    s = find_file(
        searchstr,
        exp_type=exp_type,
        expno=nodename,
        postproc=postproc,
        lookup=lookup_table,
        fl=fl,
    )
    s.ift("t2")
    # {{{ Different versions of the pulse program (notably bruker vs. spincore)
    #    phase cycle the pulses with different numbers of steps, so we need to deal
    #    with that here
    if "ph2" in s.dimlabels:
        s = s["ph2", -2]["ph1", 1]
    else:
        s = s[
            "ph1", 1
        ]  # the spincore version only has one phasing dimension 'nPhaseSteps'
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    # }}}
    # {{{ centering CPMG echo
    center = hermitian_function_test(s)
    s = center_echo(s, center, fl=fl)
    logger.debug(strm(ndshape(s)))
    fl.next("centered echo")
    fl.image(s.C.setaxis("tE", "#").set_units("tE", "scan #"))
    # }}}
    # {{{select echo decay fit function
    s.ft("t2")
    fl.next("selected coherence")
    fl.image(s.C.setaxis("tE", "#").set_units("tE", "scan #"))
    s = s["t2":f_range]
    s = s.C.sum("t2")
    # {{{ (FROM REVIEW) what is going on here???
    if "indirect" in s.dimlabels:
        CPMG = s["indirect", -1]
    else:
        CPMG = s
    # }}}
    fl.next("decay curve")
    fl.plot(CPMG, "o")
    fit_CPMG = fitdata(CPMG)
    M0, R2, vd = sp.symbols("M_0 R_2 tE", real=True)
    fit_CPMG.functional_form = (M0) * sp.exp(-vd * R2)
    logger.info(strm("Functional Form", fit_CPMG.functional_form))
    logger.info(strm("Functional Form", fit_CPMG.functional_form))
    fit_CPMG.fit()
    logger.info(strm("output:", fit_CPMG.output()))
    logger.info(strm("latex:", fit_CPMG.latex()))
    T2 = 1.0 / fit_CPMG.output("R_2")
    fl.next("CPMG fit of T1CPMG for free D2O")
    fl.plot(CPMG, "o", label="data")
    fl.plot(fit_CPMG.eval(100), label="fit")
    logger.info(strm("T2 IS:", T2))
    # }}}
    # {{{saving figure
    save_fig = False
    if save_fig:
        savefig(
            "20200108_CPMG_trials.png",
            transparent=True,
            bbox_inches="tight",
            pad_inches=0,
            legend=True,
        )
    fl.show()
    # }}}
