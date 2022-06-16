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
excluded_pathways = [(0, 0), (0,2), (0, 3)]
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
    opt_shift,sigma, my_mask = correl_align(data, indirect_dim = 'nScans',
            signal_pathway=signal_pathway, sigma = 1500)
    data.ift('t2')
    data *= np.exp(-1j*2*pi*opt_shift*data.fromaxis('t2'))
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
    #fl.show();quit()
    # }}}
    # {{{Normalization
    int_frq_slice_start = integrate_limits(select_pathway(data['nScans',0].C, signal_pathway),cutoff = 0.1)
    int_frq_slice_end = integrate_limits(select_pathway(data['nScans',-1].C,signal_pathway),cutoff = 0.1)
    int_slice = (int_frq_slice_start[0],int_frq_slice_end[-1])
    spacing = int_slice[-1]-int_slice[0]
    offres_start = int_slice[-1]+500
    offres_end = offres_start + (int_slice[-1]-int_slice[0])
    offres_slice = (offres_start,offres_end)
    s_integral = select_pathway(data['t2':int_slice].C, signal_pathway).integrate('t2')
    avg_d = s_integral.C.mean().real.item()
    s_integral /= avg_d
    data /= avg_d
    # }}}
    integral_diff_sq=abs(s_integral.C.real - s_integral.C.mean().real)**2
    f= data.getaxis('t2')
    df = f[1]-f[0]
    fl.next('limits')
    fl.plot(select_pathway(data,signal_pathway))
    plt.axvline(int_slice[0],color='blue')
    plt.axvline(int_slice[-1],color='blue')
    plt.axvline(offres_slice[0],color='red')
    plt.axvline(offres_slice[-1],color='red')
    s_intregion = select_pathway(data['t2':int_slice].C,signal_pathway)
    s_offres = select_pathway(data['t2':offres_slice].C,signal_pathway)
    s_inactive = data['t2':int_slice].C
    s_inactive['ph1',1] = 0
    N = ndshape(s_intregion.C)['t2']
    #ndshape of the t2 of the off resonance is consistently 
    #1 less than the on resonance slice
    assert (ndshape(s_offres)['t2']+1) == N
    inactive_var = s_inactive.run(var,'t2')
    s_inactive /= 2 # equiv to variance of real datapoints
    s_offres = s_offres.real # only keep the real part
    s_offres.run(var,'t2')
    # {{{ convert from error of datapoints to propagated
    #     error for integral
    s_inactive *= df **2
    s_inactive *= N
    s_offres *= df **2
    s_offres *= N
    s_offres_avg = s_offres.C.mean().item()
    # }}}
    s_intregion.integrate('t2')
    s_intregion.run(var,'nScans')
    fl.next('variances')
    plt.axhline(y = s_intregion.data, color='k',label='variance of integral')
    fl.plot(s_offres,'o',color='red',label = 'variance from off res slice')
    plt.axhline(y = s_offres_avg,color='red', label = 'averaged off res')
    for j in range(ndshape(s_inactive)['ph1']):
        if j == 1: #this is the active pathway so I skip
            pass
        else:
            fl.plot(s_inactive['ph1',j],'o', label='variance of ph1 = %d'%j)
    avg_inactive = 1/(ndshape(s_inactive)['ph1']-1)*s_inactive.sum('ph1').mean().item()
    plt.axhline(y = avg_inactive,color = 'green',label='averaged inactive')
    fl.plot(integral_diff_sq,'o',color='k',label = 'integral_diff_sq')
    print(avg_inactive)
    fl.show();quit()
    print("variance of integral",s_intregion)
    print("variance of integral, predicted based on propagation of offres",s_offres)
    print("variance of integral, predicted based on propagation of inactive:",s_inactive)
    for j in range(ndshape(s_inactive)['ph1']):
        print(s_inactive.getaxis('ph1')[j],'-->',s_inactive['ph1',j])
    # in the following the denominator doesn't include the
    # active pathway
    print('and averaged:',
            1/(ndshape(s_inactive)['ph1']-1)*s_inactive.sum('ph1'))
