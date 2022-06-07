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
signal_pathway = {"ph1": 1}
f_range = (-500, 1.5e3)
excluded_pathways = [(0, 0), (0, 3)]
colors = ["r", "darkorange", "gold", "g", "c", "b", "m", "lightcoral"]
for thisfile, exp_type, nodename in [
    ("220606_150mM_TEMPOL_200scan",
        "ODNP_NMR_comp/Echoes", "echo")
    ]:
    # {{{processing data
    data = find_file(
        thisfile,
        exp_type=exp_type,
        expno=nodename,
        postproc="spincore_Hahn_echoph_v2",
        lookup=lookup_table,
    )
    data *= -1 #frequency drift causes inverted data
    data.ift("t2")
    # {{{DC offset correction
    data.ift(["ph1"])
    t_rx = data.C.getaxis('t2')[-1]
    t_rx = data["t2":(t_rx*0.75, None)]
    t_rx = t_rx.mean(['t2'])
    data -= t_rx
    data.ft(["ph1"])
    # }}}
    data.ft("t2")
    fl.next('raw')
    fl.image(data)
    data = data["t2":f_range]
    # {{{Phase corrections
    data.ift("t2")
    best_shift = hermitian_function_test(select_pathway(data.C.mean('nScans'), signal_pathway))
    data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    data /= zeroth_order_ph(select_pathway(data,signal_pathway))
    data.ft("t2")
    # }}}
    data.ift("t2")
    data = data["t2" : (0,None)]
    data["t2":0] *= 0.5
    data.ft("t2")
    data.reorder(["ph1", "nScans", "t2"])
    # }}}
    # {{{Normalization
    frq_slice = integrate_limits(select_pathway(data.C, signal_pathway))
    s_integral = select_pathway(data['t2':frq_slice].C, signal_pathway).integrate('t2')
    avg_d = s_integral.C.mean().real.item()
    s_integral /= avg_d
    data /= avg_d
    # }}}
    error_pathway = (
        set(
            (j for j in range(ndshape(data)["ph1"]))
        )
        - set(excluded_pathways)
        - set([(signal_pathway["ph1"])])
    )
    error_pathway = [{"ph1": j} for j in error_pathway]
    # {{{Making lists for all individual inactive pathways to get error
    # associated with each one
    var_lst = []
    std_lst = []
    avg_var_lst = []
    avg_std_lst = []
    for thispathway in error_pathway:
        s_thisint, frq_slice_check = integral_w_errors(
            data.C,
            signal_pathway,
            [thispathway],
            indirect="nScans",
            return_frq_slice=True,
            fl = fl
        )
        assert all(frq_slice_check == frq_slice)
        std = s_thisint.get_error()
        var = std**2
        var_lst.append(var)
        std_lst.append(std)
        avg_var = var.mean().item()
        avg_std = sqrt(avg_var)
        avg_var_lst.append(avg_var)
        avg_std_lst.append(avg_std)
    # }}}
# {{{ Calculating propagated error averaged over all inactive CTs (as the
    #     function is meant to be called)
    data_int, frq_slice = integral_w_errors(
        data.C,
        signal_pathway,
        error_pathway,
        indirect="nScans",
        return_frq_slice=True,
    )
    inactive_std = data_int.get_error()
    inactive_var = inactive_std**2
    avg_inactive_var = inactive_var.mean().item()
    avg_inactive_std = sqrt(avg_inactive_var)
    # }}}
    # {{{ Calculating propagated error along active CT on noise slice
    active_int = active_propagation(data, signal_pathway, indirect="nScans")
    active_std = active_int.get_error()
    active_var = active_std**2
    avg_active_var = active_var.mean().item()
    avg_active_std = sqrt(avg_active_var)
    # }}}
    # {{{ Calculating the std dev -- error associated with the integrals
    expected = s_integral.run(real).mean('nScans', std=True)
    expected_std = expected.get_error()
    expected_var = expected_std**2
    # }}}
    # {{{ Plotting Errors
    fl.next("Comparison of propagated error - Real Data", legend=True)
    for i in range(len(var_lst)):
        fl.plot(
            var_lst[i],
            "o",
            color=colors[i],
            label="on excluded path of %s" % error_pathway[i],
        )
    fl.plot(active_var, "x", alpha=0.5,
            label="from integration of\nactive CT in noise slice")
    fl.plot(inactive_var, "o", color="brown", alpha=0.5,
        label="from all inactive CTs")
    axhline(
        y=avg_active_var,
        linestyle="--",
        label="averaged variance\nfrom active CT in noise slice",
    )
    axhline(
        y=avg_inactive_var,
        linestyle="--",
        color="brown",
        label="averaged variance\nfrom all inactive CTs",
    )
    axhline(
        y=expected_var,
        c="k",
        linestyle="-",
        label="expected variance",
    )
    ylabel('variance Associated')
    xlabel('nScans')
    fl.next('comparison of standard deviations',legend=True)
    for i in range(len(std_lst)):
        fl.plot(
            std_lst[i],
            "o",
            color=colors[i],
            label="on excluded path of %s" % error_pathway[i],
        )
    fl.plot(active_std, "x", 
            label="propagated error from integration of\nactive CT in noise slice")
    fl.plot(inactive_std, "o", color="brown",
        label="Error from all inactive CTs")
    for i in range(len(std_lst)):
        axhline(
            y=avg_std_lst[i],
            linestyle=":",
            color=colors[i],
            label="averaged %s" % error_pathway[i],
        )
    axhline(
        y=avg_active_std,
        linestyle="--",
        label="averaged propagated error\nfrom active CT in noise slice",
    )
    axhline(
        y=avg_inactive_std,
        linestyle="--",
        color="brown",
        label="averaged propagated error\nfrom all inactive CTs",
    )
    axhline(
        y=expected_std,
        c="k",
        linestyle="-",
        label="traditional std using numpy",
    )
    xlabel('nScans')
    ylabel('Associated Std')
    #}}}
    print("**** *** ***")
    print("estimated noise of the integral:",expected_std)
    print("**** *** ***")
    print("**** *** ***")
    print("propagated error along inactive pathways on noise slice (integral_w_errors):",avg_inactive_std)
    print("**** *** ***")
    print("**** *** ***")
    print("propagated error along active pathways on noise slice (active_propagation):",avg_active_std)
    print("**** *** ***")
    plt.axis("tight")
    ax = plt.gca()
    lims = list(ax.get_ylim())
    lims[0] = 0
    ax.set_ylim(lims)
    plt.legend()
    fl.show()

