from pyspecdata import *
from pylab import *
from matplotlib import *
from pyspecProcScripts import *
from pyspecProcScripts.correlation_alignment import correl_align
import numpy as np
fl = figlist_var()
signal_pathway = {'ph1': 1, 'ph2':0}
excluded_pathways = [(0,0),(0,3)]
for thisfile,exp_type,nodename in [
        ('201113_TEMPOL_capillary_probe_16Scans_noModCoil','ODNP_NMR_comp/test_equipment','signal')
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename)
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = 8
    SW_kHz =s.get_prop('acq_params')['SW_kHz']
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
    t_range=(0,0.05)
    f_range = (-0.11e3,0.12e3)
    s.ift(['ph1','ph2'])
    t_rx = (t_range[-1]/4)*3
    rx_offset_corr = s['t2':(t_rx,None)]
    rx_offset_corr = rx_offset_corr.data.mean()
    s -= rx_offset_corr
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
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',1]['ph2',0] = 1
    phasing.ift(['ph1','ph2'])
    s /= phasing
    fl.next('ph, freq -- apply reciever phase')
    fl.image(s)
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
    s /= ph0
    s.ft(['ph1','ph2'])
    fl.next('after zeroth order phasing applied')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    #{{{alignment
    #s.ift(['ph1','ph2'])
    #s.ft('t2')
    #s_aligned,opt_shift,sigma = correl_align(s,indirect_dim='nScans',
    #        ph1_selection = signal_pathway['ph1'],
    #        ph2_selection = signal_pathway['ph2'],sigma=50)
    #s = s_aligned
    #s.ift('t2')
    #s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    #s.ft('t2')
    #s = s_aligned
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
    s.ift('t2')
    s = s['t2':(0,t_range[-1])]
    s['t2':0] *= 0.5
    s.ft('t2')
    fl.next('FID sliced')
    fl.image(s)
    s.reorder(['ph1','ph2','nScans','t2'])
    s /= s.data.mean()
    data = s.C  
    error_pathway = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    s_int,frq_slice = integral_w_errors(s,signal_pathway,error_pathway,
            indirect='nScans',fl=fl,return_frq_slice=True)
    x = s_int.get_error()
    x[:] /= 2
    fl.next('comparison of std')
    avg_s_int = s_int.get_error().mean().item()
    s_int1 = avg_s_int
    s_int1 /= s_int.data.argmax().item()
    fl.plot(s_int.get_error(),'o',label='returned error from integral_w_errors')
    axhline(y=avg_s_int,c='red',linestyle=":",label='averaged returned error from integral_w_errors')
    data1 =data.C
    data = select_pathway(data,signal_pathway)
    data.integrate('t2')
    data1 = data1['t2':frq_slice]
    data2 = data1.C
    data1=select_pathway(data1,signal_pathway)
    data2=select_pathway(data2,{'ph1':0,'ph2':0})
    data2 /= data2.argmax('t2').item()
    data1 /= data1.argmax('t2').item()

    data2.integrate('t2')
    data1.integrate('t2')
    after_slice_error = data1.real.run(np.std,'nScans')
    off_CT_error = data2.real.run(np.std,'nScans')
    axhline(y=float(after_slice_error.data),
            c='k',linestyle=":",label='error associated with CT pathway of signal')
    axhline(y=float(off_CT_error.data),
            c='blue',linestyle=":",label='error associated with inactive CT pathway')
    print('error from integral w errors',avg_s_int)
    print('error associated with CT pathway of signal', float(after_slice_error.data))
    print('error associated with inactive CT pathway',float(off_CT_error.data))
    plt.axis('tight')
    ax = plt.gca()
    lims = list(ax.get_ylim())
    lims[0] = 0
    ax.set_ylim(lims)
    plt.legend()
    print("error from integral w errors",avg_s_int)
    print("error associated with CT of signal", float(after_slice_error.data))
    print("error associated with inactive CT pathway",float(off_CT_error.data))
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

