""" Validate Inactive CT Error Using Real Data
==============================================

Estimates the error of the integral of an actual data set of a standard echo
experiment. Three methods of acquiring the error associated with the data are 
compared:

    -   Taking an area along the active coherence transfer (CT) pathway outside of the bandwidth of the signal
        signal and propagating that error to estimate the error associated with the integral.
        (The traditional method of acquiring the error associated with a data set.)
    -   Taking the integral in the inactive CT pathways and propagating to get the error 
        associated with the integral in the active CT.
    -   Taking the standard deviation of many integrals determined by
        integrating over the signal bandwidth of the active CT pathway.
        (Best method when many scans are available)
    
Demonstrates that by propagating the error of the integral in the inactive CTs we still
get a reasonable error within the limits of traditional methods.
"""
from pyspecdata import *
from pylab import *
from matplotlib import *
from pyspecProcScripts import *
from pyspecProcScripts.correlation_alignment import correl_align
import numpy as np
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 4

fl = figlist_var()
signal_pathway = {"ph1": 1, "ph2": 0}
t_range = (0, 0.05) # must start at 0 for an FID
f_range = (-1e3, 1e3)
excluded_pathways = [(0, 0), (0, 3)]
colors = ["r", "darkorange", "gold", "g", "c", "b", "m", "lightcoral"]
for thisfile, exp_type, nodename in [
    ("201113_TEMPOL_capillary_probe_16Scans_ModCoil",
        "ODNP_NMR_comp/old", "signal")
    ]:
    # {{{processing data
    s = find_file(
        thisfile,
        exp_type=exp_type,
        expno=nodename,
        postproc="spincore_echo_v1",
        lookup=lookup_table,
    )
    s.ift("t2")
    # {{{DC offset correction
    s.ift(["ph1", "ph2"])
    t_rx = s.C.getaxis('t2')[-1]
    t_rx = s["t2":(t_rx*0.75, None)]
    t_rx = t_rx.mean(['t2'])
    s -= t_rx
    s.ft(["ph1", "ph2"])
    # }}}
    s.ft("t2")
    s = s["t2":f_range]
    # {{{Phase corrections
    s.ift("t2")
    best_shift = hermitian_function_test(select_pathway(s.C.mean('nScans'), signal_pathway))
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    s /= zeroth_order_ph(select_pathway(s['t2':0],signal_pathway))
    s.ft("t2")
    # }}}
    s.ift("t2")
    s = s["t2" : (0,None)]
    s["t2":0] *= 0.5
    s.ft("t2")
    s.reorder(["ph1", "ph2", "nScans", "t2"])
    # }}}
    # {{{Normalization
    frq_slice = integrate_limits(select_pathway(s.C, signal_pathway))
    s_integral = select_pathway(s['t2':frq_slice].C, signal_pathway).integrate('t2')
    avg_d = s_integral.C.mean().real.item()
    s_integral /= avg_d
    s /= avg_d
    # }}}
    error_pathway = (
        set(
            ((j, k) for j in range(ndshape(s)["ph1"]) for k in range(ndshape(s)["ph2"]))
        )
        - set(excluded_pathways)
        - set([(signal_pathway["ph1"], signal_pathway["ph2"])])
    )
    error_pathway = [{"ph1": j, "ph2": k} for j, k in error_pathway]
    # {{{Making lists for all individual inactive pathways to get error
    # associated with each one
    error_lst = []
    avg_error_lst = []
    for thispathway in error_pathway:
        s_thisint, frq_slice_check = integral_w_errors(
            s,
            signal_pathway,
            [thispathway],
            indirect="nScans",
            return_frq_slice=True,
        )
        assert all(frq_slice_check == frq_slice)
        error = s_thisint.get_error()
        avg_error = error.mean().item()
        error_lst.append(error)
        avg_error_lst.append(avg_error)
    # }}}
    # {{{ Calculating propagated error averaged over all inactive CTs (as the
    #     function is meant to be called)
    averaged_inactive, frq_slice = integral_w_errors(
        s.C,
        signal_pathway,
        error_pathway,
        indirect="nScans",
        return_frq_slice=True,
    )
    averaged_inactive_error = averaged_inactive.get_error()
    avg_avg_error = averaged_inactive_error.mean().item()
    # }}}
    # {{{ Calculating propagated error along active CT on noise slice
    active_error = active_propagation(s, signal_pathway, offset = 700, indirect="nScans")
    avg_active_error = active_error.get_error()
    avg_avg_active_err = avg_active_error.mean().item()
    # }}}
    # {{{ Calculating the std dev -- error associated with the integrals
    expected = s_integral.run(real).mean('nScans', std=True)
    s_integral = expected.get_error()
    # }}}
    # {{{ Plotting Errors
    fl.next("comparison of propagated error", legend=True)
    for i in range(len(error_lst)):
        fl.plot(error_lst[i], "o", color=colors[i], alpha=0.5, 
            label="Error on excluded path of \n %s" % error_pathway[i])
    fl.plot(avg_active_error, "x", alpha=0.5,
            label="propagated error from active CT\nin noise slice")
    fl.plot(
        averaged_inactive_error, "o", color="brown", alpha=0.5,
        label="averaged propagated error\nfrom all inactive CTs")
    for i in range(len(error_lst)):
        axhline(y=avg_error_lst[i], linestyle=":", color=colors[i],
            label="averaged propagated error over \n %s" % error_pathway[i])
    axhline(
        y=avg_avg_active_err,
        linestyle="--",
        label="averaged propagated error\nfrom active CT in noise slice",
    )
    axhline(
        y=avg_avg_error,
        linestyle="--",
        color="brown",
        label="averaged propagated error\nfrom all inactive CTs",
    )
    axhline(
        y=s_integral,
        c="k",
        linestyle="-",
        label="traditional std using numpy",
    )
    print("**** *** ***")
    print("estimated noise of the integral:",s_integral)
    print("**** *** ***")
    print("**** *** ***")
    print("propagated error along inactive pathways on noise slice (integral_w_errors):",avg_avg_error)
    print("**** *** ***")
    print("**** *** ***")
    print("propagated error along active pathways on noise slice (active_propagation):",avg_avg_active_err)
    print("**** *** ***")
    # }}}
    plt.axis("tight")
    ax = plt.gca()
    plt.legend()
    fl.show()

