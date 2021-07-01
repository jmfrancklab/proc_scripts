from pyspecdata import *
from pylab import *
from matplotlib import *
from pyspecProcScripts import *
from pyspecProcScripts.correlation_alignment import correl_align
import numpy as np
fl = figlist_var()
signal_pathway = {'ph1': 1, 'ph2':0}
t_range=(0,0.05)
f_range = (-0.11e3,0.12e3)
excluded_pathways = [(0,0),(0,3)]
colors = ['r','darkorange','gold','g','c','b','m','lightcoral']
for thisfile,exp_type,nodename in [
        ('201113_TEMPOL_capillary_probe_16Scans_noModCoil','ODNP_NMR_comp/test_equipment','signal')
        ]:
#{{{processing data    
    s = find_file(thisfile,exp_type=exp_type,expno=nodename)
    nScans = s.get_prop('acq_params')['nScans']
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph1','ph2'])
    s.setaxis('nScans',r_[0:nScans])
    s.set_units('t2','s')
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    fl.next('raw data')
    fl.image(s)
    s.ift('t2')
    fl.next('raw data time domain')
    fl.image(s)
    #{{{DC offset correction
    s.ift(['ph1','ph2'])
    t_rx = (t_range[-1]/4)*3
    s -= s['t2':(t_rx,None)].data.mean()  
    s.ft('t2')
    s.ft(['ph1','ph2'])
    #}}}
    s = s['t2':f_range]
    fl.next('freq domain')
    fl.image(s)
    #{{{Phase corrections
    s.ift('t2')
    best_shift,window_size = hermitian_function_test(select_pathway(s,signal_pathway))
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0})
    fl.next('After hermitian phase correction')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    ph0 = select_pathway(s,signal_pathway)['t2':0]
    ph0 /= abs(ph0)
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
    s /= ph0
    fl.next('after zeroth order phasing applied')
    s.ft('t2')
    fl.image(s)
    #}}}
    s.ift('t2')
    #{{{alignment
    #s.ift(['ph1','ph2'])
    #s.ft('t2')
    #opt_shift,sigma = correl_align(s,indirect_dim='nScans',
    #        ph1_selection = signal_pathway['ph1'],
    #        ph2_selection = signal_pathway['ph2'],sigma=50)
    #s.ift('t2')
    #s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    #s.ft('t2')
    #fl.basename=None
    #fl.next(r'after correlation, $\varphi$ domain')
    #fl.image(s)
    #s.ift('t2')
    #s.ft(['ph1','ph2'])
    #fl.next('after correl - time domain')
    #fl.image(s)
    #s.ft('t2')
    #fl.next('after correl - freq domain')
    #fl.image(s)
    ##}}}
    #s.ift('t2')
    s = s['t2':(0,t_range[-1])]
    s['t2':0] *= 0.5
    s.ft('t2')
    fl.next('FID sliced')
    fl.image(s)
    s.reorder(['ph1','ph2','nScans','t2'])
    #}}}
    
    #{{{Normalization
    frq_slice = integrate_limits(select_pathway(s,signal_pathway),fl=fl)
    d = s.C['t2':frq_slice]
    d = select_pathway(d,signal_pathway)
    d.integrate('t2')
    avg_d = d.C.mean().item()
    d /= avg_d
    s /= avg_d
    #}}}
    
    #{{{integral w errors
    error_pathway = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    #{{{Making lists for all individual inactive pathways to get error
    # associated with each one
    s_int_lst = []
    ph_lst = [[{'ph1':0,'ph2':0}],
        [{'ph1':0,'ph2':1}],
        [{'ph1':1,'ph2':1}],
        [{'ph1':2,'ph2':0}],
        [{'ph1':2,'ph2':1}],
        [{'ph1':3,'ph2':0}],
        [{'ph1':3,'ph2':1}]]
    error_lst = []
    avg_error_lst = []
    for i in range(len(ph_lst)):
        test = s.get_ft_prop('t2')
        if s.get_ft_prop('t2') is False:
            s.ft('t2')
        s_int,frq_slice = integral_w_errors(s, signal_pathway, ph_lst[i],
                    indirect='nScans', fl=fl, return_frq_slice=True)
        error = s_int.get_error()
        error[:] /= 2
        avg_error = error.mean().item()
        s_int_lst.append(s_int)
        error_lst.append(error)
        avg_error_lst.append(avg_error)
    #}}}

    #{{{Calculating propagated error averaged over all inactive CTs
    averaged_inactive, frq_slice = integral_w_errors(s, signal_pathway,
            error_pathway, indirect='nScans', fl=fl, return_frq_slice=True)
    averaged_inactive_error = averaged_inactive.get_error()
    averaged_inactive_error[:] /= 2
    avg_avg_error = averaged_inactive_error.mean().item()
    #}}}

    #{{{Calculating propagated error along active CT on noise slice
    active_error,N,df = active_propagation(s, signal_pathway, indirect='nScans',fl=fl)
    active_error[:] /= 2
    avg_active_error = active_error.mean().item()
    #}}}
    #{{{Calculating the std dev -- error associated with the integrals
    numpy_s_int = d.run(np.std,'nScans')
   #}}} 
    #{{{Plotting Errors
    fl.next('comparison of std')
    for i in range(len(s_int_lst)):
        fl.plot(error_lst[i],'o',color=colors[i],label = 'on excluded path of %s'%ph_lst[i])
    fl.plot(active_error,'x',
            label='propagated error from active CT in noise slice')
    fl.plot(averaged_inactive_error,'o',color='brown',label='averaged propagated error from all inactive CTs')
    for i in range(len(s_int_lst)):
        axhline(y=avg_error_lst[i], linestyle=":", color=colors[i],
                label = "averaged %s"%ph_lst[i])
    axhline(y=avg_active_error,linestyle=":", label='averaged propagated error from active CT in noise slice')
    axhline(y=avg_avg_error,linestyle=":",color='brown',label='averaged average propagated error from inactive CTs')
    axhline(y=numpy_s_int.data,c='k',
            linestyle=":",label='std dev - error associated with the integrals')
    #}}}
    plt.axis('tight')
    ax = plt.gca()
    lims = list(ax.get_ylim())
    lims[0] = 0
    ax.set_ylim(lims)
    plt.legend()
    fl.show();quit()

