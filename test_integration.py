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
for thisfile,exp_type,nodename in [
        ('201113_TEMPOL_capillary_probe_16Scans_noModCoil','ODNP_NMR_comp/test_equipment','signal')
        ]:
#{{{processing data    
    s = find_file(thisfile,exp_type=exp_type,expno=nodename)
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph1','ph2'])
    s.setaxis('nScans','#')
    s.set_units('t2','s')
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    fl.next('raw data')
    fl.image(s)
    s.ift('t2')
    fl.next('raw data time domain')
    fl.image(s)
    s.ift(['ph1','ph2'])
    t_rx = (t_range[-1]/4)*3
    s -= s['t2':(t_rx,None)].data.mean()  # DC offset correction
    s.ft('t2')
    s.ft(['ph1','ph2'])
    s = s['t2':f_range]
    fl.next('freq domain')
    fl.image(s)
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
    #fl.next('FID sliced')
    #fl.image(s)
    s.reorder(['ph1','ph2','nScans','t2'])
    #}}}
    d = s.C.integrate('t2')
    s /= d.data.mean()

    #{{{integral w errors
    error_pathway = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    s_int_0, active_error, frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':0,'ph2':0}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_1, active_error,frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':0,'ph2':1}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_2, active_error,frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':1,'ph2':1}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_3, active_error,frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':2,'ph2':0}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_4, active_error,frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':2,'ph2':1}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_5, active_error,frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':3,'ph2':0}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_6, active_error,frq_slice = integral_w_errors(s,signal_pathway,[{'ph1':3,'ph2':1}],
            indirect='nScans',fl=fl,return_frq_slice=True)
    s_int_all, active_error,frq_slice = integral_w_errors(s,signal_pathway,error_pathway,
            indirect='nScans',fl=fl,return_frq_slice=True)
    S = s_int_all.get_error()
    #np.std(s_int.data)
    S[:] /= 2
    x = s_int_0.get_error()
    x[:] /= 2
    y = active_error.get_error()
    y[:] /= 2
    z = s_int_1.get_error()
    z[:] /= 2
    a = s_int_2.get_error()
    a[:] /= 2
    b = s_int_3.get_error()
    b[:] /= 2
    c = s_int_4.get_error()
    c[:] /= 2
    d = s_int_5.get_error()
    d[:] /= 2
    e = s_int_6.get_error()
    e[:] /= 2
    fl.next('comparison of std')
    avg_s_int_0 = s_int_0.get_error().mean().item()
    avg_s_int_1 = s_int_1.get_error().mean().item()
    avg_s_int_2 = s_int_2.get_error().mean().item()
    avg_s_int_3 = s_int_3.get_error().mean().item()
    avg_s_int_4 = s_int_4.get_error().mean().item()
    avg_s_int_5 = s_int_5.get_error().mean().item()
    avg_s_int_6 = s_int_6.get_error().mean().item()
    avg_s_all = s_int_all.get_error().mean().item()
    avg_s_int_active = active_error.get_error().mean().item()
    fl.plot(s_int_0.get_error(),'o',label='propagated error for ph1:0, ph2:0')
    fl.plot(s_int_1.get_error(),'o',label='propagated error for ph1:0, ph2:1')
    fl.plot(s_int_2.get_error(),'o',label='propagated error for ph1:1, ph2:1')
    fl.plot(s_int_3.get_error(),'o',label='propagated error for ph1:2, ph2:0')
    fl.plot(s_int_4.get_error(),'o',label='propagated error for ph1:2, ph2:1')
    fl.plot(s_int_5.get_error(),'o',label='propagated error for ph1:3, ph2:0')
    fl.plot(s_int_6.get_error(),'o',label='propagated error for ph1:3, ph2:1')
    fl.plot(s_int_all.get_error(),'x',label='propagated error over all excluded paths')
    fl.plot(active_error.get_error(),'o',label='propagated error from active CT in noise slice')
    axhline(y=avg_s_int_0,c='red',linestyle=":",label='averaged propagated error for ph1:0,ph2:0')
    axhline(y=avg_s_int_1,c='orange',linestyle=":",label='averaged propagated error for ph1:0,ph2:1')
    axhline(y=avg_s_int_2,c='yellow',linestyle=":",label='averaged propagated error for ph1:1,ph2:1')
    axhline(y=avg_s_int_3,c='green',linestyle=":",label='averaged propagated error for ph1:2,ph2:0')
    axhline(y=avg_s_int_4,c='cyan',linestyle=":",label='averaged propagated error for ph1:2,ph2:1')
    axhline(y=avg_s_int_5,c='blue',linestyle=":",label='averaged propagated error for ph1:3,ph2:0')
    axhline(y=avg_s_int_6,c='blue',linestyle=":",label='averaged propagated error for ph1:3,ph2:1')
    axhline(y=avg_s_all,c='magenta',linestyle=":",label='averaged propagated error for all excluded paths')

    axhline(y=avg_s_int_active,c='blue',linestyle=":",label='averaged propagated error from active CT in noise slice')
    
    #}}}
    #{{{numpy stds
    np_std = s_int_all.real.run(np.std('nScans'))
    axhline(y=float(np_std.data),c='k',
            linestyle=":",label='std noise slice of CT')
    #}}}
    plt.axis('tight')
    ax = plt.gca()
    lims = list(ax.get_ylim())
    lims[0] = 0
    ax.set_ylim(lims)
    plt.legend()
    fl.show();quit()
    fl.next('diagnostic 1D plot')
    fl.plot(s['nScans',:]['ph1',signal_pathway['ph1']]['ph2',signal_pathway['ph2']].real,alpha=0.4)
    axvline(x=frq_slice[0],c='k',linestyle=":",alpha=0.8)
    axvline(x=frq_slice[-1],c='k',linestyle=":",alpha=0.8)
    fl.next('integrated for error')
    fl.plot(s_int,'.',capsize=6,label='integral with error')
    #data = data['t2':frq_slice]
    data = select_pathway(data,signal_pathway)
    data.integrate('t2')

    new_error = data.real.run(np.std,'nScans')
    #new_error = data.real.C.mean('nScans', std=True).get_error()
    fl.plot(data,'o',label='data')
    #print("new_error, data",ndshape(new_error), ndshape(data))
    data.set_error(new_error.data)
    fl.plot(data,'.',capsize=6,label='numpy std')
    fl.show();quit()

