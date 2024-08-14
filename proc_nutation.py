from pylab import *
from pyspecdata import *
from scipy.optimize import minimize
from pyspecProcScripts import (
    hermitian_function_test,
    zeroth_order_ph,
    lookup_table,
    fl_mod,
    DCCT,
)
from sympy import symbols
from numpy import *

fl = fl_mod()
t2 = symbols("t2")
logger = init_logging("info")
max_kHz = 200
for searchstr, exp_type, nodename, postproc, freq_slice, t_slice in [
        # PR: either choose to process all, or only one
    # ['201211_Ni_sol_probe_nutation_1','nutation','nutation',
    #    'spincore_nutation_v1',(-5000,13000)],
    # ['210302_210302_Ni_cap_probe_nutation_1','nutation','nutation',
    #    'spincore_nutation_v1',(-4e3,1.2e4)]
    [
        "210409_Ni_sol_probe_nutation_1",
        "nutation",
        "nutation",
        "spincore_nutation_v3",
        (11e3, 40e3),
        (0.5e-3, 1.5e-3),
    ]
]:
    s = find_file(
        searchstr,
        exp_type=exp_type,
        expno=nodename,
        postproc=postproc,# PR why do you need to manually specify the postproc? that should not be needed!
        lookup=lookup_table,
        fl=fl,
    )
    s = s["t2":freq_slice]
    fl.next("data")
    fl.image(s)
    # {{{ PR why are we doing the following -- is it needed?? If it is, comment on this, and also include a slide in your pptx
    s.ift("t2")
    rx_offset_corr = s[
        "t2":(0.75e-3, None)
    ]  # should be quarter of t_slice onward
    rx_offset_corr = rx_offset_corr.mean(["t2"])
    s -= rx_offset_corr
    s.ft("t2")
    fl.next("After rx offset correction")
    fl.image(s)
    # }}}
    if "amp" in s.dimlabels:
        plen = s.get_prop("acq_params")["p90_us"] * 1e-6
        logger.info(strm("pulse length is:", plen))
    s.ift("t2")
    # {{{ PR why would you not just use the fid_from_echo here??
    # {{{ do the centering before anything else!
    # in particular -- if you don't do this before convolution, the
    # convolution doesn't work properly!
    print(ndshape(s))
    best_shift, window = hermitian_function_test(s["ph1", 1]["ph2", -2])
    logger.info("best shift is:", best_shift)
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    fl.next("with time shift")
    fl.image(s)
    fl.next("freq domain after time correction")
    s.ft("t2")
    fl.image(s)
    s.ift("t2")  # make sure everything is in the same domain
    # }}}
    # }}}
    fl.next("t domain centered")
    fl.image(s)
    # {{{zeroth order phasing
    s.ift(["ph1", "ph2"])
    phasing = s["t2", 0].C
    phasing.data *= 0
    phasing.ft(["ph1", "ph2"])
    phasing["ph1", 1]["ph2", -2] = 1
    phasing.ift(["ph1", "ph2"])
    fl.next("zeroth order corrected")
    ph0 = s["t2":0] / phasing
    ph0 /= abs(ph0)
    s /= ph0
    s.ft(["ph1", "ph2"])
    fl.next("phased")
    ph0 = zeroth_order_ph(s["t2":0], fl=None)
    s /= ph0
    s.ft("t2")
    fl.image(s)
    # }}}
    # {{{ selecting coherence and convolving
    s = s["ph1", 1]["ph2", -2]
    fl.next("select $\\Delta p$")
    fl.image(s)
    # }}}
    if "amp" in s.dimlabels:
        s.setaxis("amp", lambda x: x * plen)
        s.set_units("amp", "s")
        ind_dim = "\\tau_p a"
        s.rename("amp", ind_dim)
    elif "p_90" in s.dimlabels:
        ind_dim = "p_90"
    else:
        raise ValueError("not sure what the indirect dimenison is!!")
    fl.image(s, human_units=False)
    fl.real_imag("phased data", s)
    fl.next("FT")
    title("FT to get $\gamma B_1/a$")
    s.ft(ind_dim, shift=True)
    fl.image(s[ind_dim : (-1e3 * max_kHz, 1e3 * max_kHz)])
    fl.next("absFT")
    title("FT to get $\gamma B_1/a$")
    fl.image(abs(s[ind_dim : (-1e3 * max_kHz, 1e3 * max_kHz)]))
    gridandtick(gca(), gridcolor=[1, 1, 1])
fl.show()
