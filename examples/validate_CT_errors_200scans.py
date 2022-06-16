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
f_range = (0, 1.1e3)
excluded_pathways = [(0,0),(0,2),(0,3)]# 0), (0, 3)]
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
    #fl.show();quit()
    #data = data["t2":f_range]
    # {{{Phase corrections
    data.ift("t2")
    best_shift = hermitian_function_test(select_pathway(data['t2':f_range].C.mean('nScans'), signal_pathway))
    print(best_shift)
    data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    data /= zeroth_order_ph(select_pathway(data['t2':0],signal_pathway))
    data.ft("t2")
    # }}}
    #{{{ Alignment??
    data.ift(['ph1'])
    #opt_shift,sigma, my_mask = correl_align(data, indirect_dim = 'nScans',
    #        signal_pathway=signal_pathway, sigma = 1500)
    data.ift('t2')
    #data *= np.exp(-1j*2*pi*opt_shift*data.fromaxis('t2'))
    data.ft(['ph1'])
    data.ft('t2')
    fl.next('aligned')
    fl.image(data)
    data.ift("t2")
    data = data["t2" : (0,None)]
    data["t2":0] *= 0.5
    data.ft("t2")
    data.reorder(["ph1", "nScans", "t2"])
    fl.next('Real Data phased - 200 scans')
    fl.image(data)
    og_data = data.C
    #fl.show();quit()
    # }}}
    #}}}
    # {{{Normalization
    int_frq_slice_start = integrate_limits(select_pathway(data['nScans',0].C, signal_pathway),cutoff = 0.1)
    int_frq_slice_end = integrate_limits(select_pathway(data['nScans',-1].C,signal_pathway),cutoff = 0.1)
    int_slice = (int_frq_slice_start[0],int_frq_slice_end[-1])
    s_integral = select_pathway(data['t2':int_slice].C, signal_pathway).integrate('t2')
    avg_d = s_integral.C.mean().real.item()
    s_integral /= avg_d
    data /= avg_d
    fl.next('Normalized Real Data')
    fl.plot(select_pathway(data.C.mean('nScans'),signal_pathway))
    # }}}
    error_pathway = (
        set(
            ((j) for j in range(ndshape(data)["ph1"]))
        )
        - set(excluded_pathways)
        - set([(signal_pathway["ph1"])])
    )
    error_pathway = [{"ph1": j} for j in error_pathway]
# {{{ Calculating propagated error averaged over all inactive CTs (as the
    #     function is meant to be called)
    fl.next('test limits last scan')
    inact_frq_slice_start = integrate_limits(
        select_pathway(data['nScans',0].C, signal_pathway),  
        cutoff=0.1, fl=fl)
    inact_frq_slice_end = integrate_limits(
            select_pathway(data['nScans',-1].C,signal_pathway),
            cutoff = 0.1)
    inact_frq_slice = (inact_frq_slice_start[0],inact_frq_slice_end[-1])
    s = data['t2':inact_frq_slice].C
    f = s.getaxis('t2')
    df = f[1] - f[0]
    errors = []
    for j in range(len(error_pathway)):
        s_forerror = select_pathway(s, error_pathway[j])
        if j == 0:
            N2 = ndshape(s_forerror)['t2']
        s_forerror -= s_forerror.C.mean_all_but(['nScans'])
        s_forerror.run(lambda x: abs(x) ** 2 / 2).mean_all_but(['nScans'])
        s_forerror *= df ** 2  # Δf
        s_forerror *= N2
    s = select_pathway(s, signal_pathway)
    inact_var = s_forerror.data
    inact_std = sqrt(s_forerror.data)
    avg_inact_var = inact_var.mean().item()
    avg_inact_std = inact_std.mean().item()
    fl.plot(select_pathway(data['t2':(-1000,2000)].C, signal_pathway))
    axvline(inact_frq_slice[-1],color = 'red')
    axvline(inact_frq_slice[0],color = 'red')
    axvline(int_slice[-1],color='blue')
    axvline(int_slice[0],color='blue')
    #fl.show();quit()
    # }}}
        # {{{Making lists for all individual inactive pathways to get error
    # associated with each one
    test = data.C
    test = test['t2':int_slice]
    test['ph1',1] = 0
    test_var = test.C.getaxis('t2').run(lambda x: abs(x) ** 2/2)
    test_var *= ((N2)*(df**2))
    print(ndshape(test_var))
    print(test_var)
    print(test_var.data)
    quit()
    # }}}

    # {{{ Calculating propagated error along active CT on noise slice
    act_frq_slice = integrate_limits(
            select_pathway(data['nScans',-1].C, signal_pathway),
            cutoff = 0.1, fl=fl)
    slice_start = act_frq_slice[-1] + 500 
    off_res_frq_slice = (slice_start,None)
    s = data['t2' : off_res_frq_slice].C
    s = s['t2',0:N2]
    s_forerror = select_pathway(s, signal_pathway)
    s_forerror -= s_forerror.C.mean_all_but(['nScans'])
    s_forerror.run(lambda x: np.real(x) ** 2).mean_all_but(['nScans'])
    s_forerror *= df ** 2
    s_forerror *= N2
    s = select_pathway(s,signal_pathway)
    active_var = s_forerror.data
    active_std = sqrt(s_forerror.data)
    avg_active_var = active_var.mean().item()
    avg_active_std = active_std.mean().item()
    axvline(s.getaxis('t2')[0],color='k')
    axvline(s.getaxis('t2')[-1],color='k')
    # }}}
    # {{{ Calculating the std dev -- error associated with the integrals
    expected = s_integral.C.run(real).mean('nScans',std=True)
    expected_std = expected.get_error()
    expected_var = expected_std**2
    integral_diff = s_integral.C.real - s_integral.C.mean().real
    integral_diff_sq = abs(s_integral.C.real - s_integral.C.mean().real)**2
    # }}}
    # {{{ Plotting Errors
    #{{{variances
    fl.next("Comparison of variances - not aligned", legend=True)
    #for i in range(len(var_lst)):
    #    fl.plot(
    #        var_lst[i],
    #        "o",
    #        color=colors[i],
    #        label="on excluded path of %s" % error_pathway[i],
    #    )
    fl.plot(active_var, "x", 
            label="propagated error from integration of\nactive CT in noise slice")
    fl.plot(inact_var, "o", color="brown",label="Error from all inactive CTs")
    #fl.plot(integral_diff_sq, 'o',color = 'k',label = 'integral diff squared')
    #for i in range(len(std_lst)):
    #    axhline(
   #         y=avg_var_lst[i],
   #         linestyle=":",
   #         color=colors[i],
   #         alpha = 0.4,
    #        label="averaged %s" % error_pathway[i],
    #    )
    axhline(y = avg_test_var,color = 'blue',label = 'averaged by hand')
    axhline(y = test_0_var,alpha = 0.5,color='red',label = 'ph1=0')
    axhline(y = test_2_var,alpha = 0.5,color='orange',label='ph1=2')
    axhline(y = test_3_var, alpha = 0.5,color='green',label = 'ph1=3')
    axhline(y=avg_active_var,linestyle="--",
        label="averaged propagated error\nfrom active CT in noise slice")
    axhline(y=avg_inact_var,linestyle="--",color="blue",
        label="averaged propagated error\nfrom all inactive CTs")
    axhline(y=expected_var, c="k",linestyle="-",
        label="traditional std using numpy")
    xlabel('nScans')
    ylabel('Associated variances')
    #}}}
    fl.next('comparison of standard deviations- not aligned',legend=True)
    #for i in range(len(std_lst)):
    #    fl.plot(
    #        std_lst[i],
    #        "o",
    #        color=colors[i],
    #        label="on excluded path of %s" % error_pathway[i],
    #    )
    fl.plot(integral_diff, 'o', color = 'k')
    fl.plot(active_std, "x", 
            label="propagated error from integration of\nactive CT in noise slice")
    fl.plot(inact_std, "o", color="brown", label="Error from all inactive CTs")
    #for i in range(len(std_lst)):
    #    axhline(
    #        y=avg_std_lst[i],
    #        linestyle=":",
    #        color=colors[i],
    #        label="averaged %s" % error_pathway[i],
    #    )
    axhline(y=avg_active_std,linestyle="--",
        label="averaged propagated error\nfrom active CT in noise slice")
    axhline(y=avg_inact_std,linestyle="--",color="orange",
        label="averaged propagated error\nfrom all inactive CTs")
    axhline(y=expected_std,c="k",linestyle="-",
        label="traditional std using numpy")
    xlabel('nScans')
    ylabel('Associated Std')
    #}}}
    print("**** *** ***")
    print("sample std:",expected_std)
    print("**** *** ***")
    print("**** *** ***")
    print("propagated error along inactive pathways on noise slice (integral_w_errors):",avg_inact_std)
    print("**** *** ***")
    print("**** *** ***")
    print("propagated error along active pathways on noise slice (active_propagation):",avg_active_std)
    print("**** *** ***")
    plt.axis("tight")
    ax = plt.gca()
    lims = list(ax.get_ylim())
    lims[0] = 0
    #ax.set_ylim(lims)
    plt.legend()
    fl.show()

