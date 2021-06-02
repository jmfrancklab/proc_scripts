from pyspecdata import *
from pylab import *
from proc_scripts import *
from proc_scripts.correlation_alignment_ODNP import correl_align
import numpy as np
fl = figlist_var()
def select_pathway(s,pathway):
    retval = s
    for k, v in pathway.items():
        retval = retval[k,v]
    return retval
signal_pathway = {'ph1': 1, 'ph2':0}
excluded_pathways = [(0,0),(0,1)]
for thisfile,exp_type,nodename in [
        ('210409_Ni_cap_probe_echo_1024','ODNP_NMR_comp/test_equipment','signal')
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
    s.chunk('nScans',['repeats','nScans'],[32,-1])
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    n_nScans = 32
    s = s.mean('repeats')
    fl.next('raw data')
    fl.image(s)
    s.ift('t2')
    fl.next('raw data time domain')
    t_range=(0,0.06)
    f_range = (-1e3,1e3)
    n_repeats=32
    all_results = ndshape(s) + (n_repeats,'repeats')
    all_results.pop("t2").pop("ph1").pop("ph2")
    all_results = all_results.alloc()
    all_results.setaxis('nScans',s.getaxis('nScans'))
    s.ift(['ph1','ph2'])
    rx_offset_corr = s['t2':(0.045,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    s.ft(['ph1','ph2'])
    s = s['t2':f_range]
    s.ift('t2')
    print("4")
    best_shift,window_size = hermitian_function_test(select_pathway(s,signal_pathway).C.convolve('t2',0.01))
    print("5")
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
    fl.next('phase corrections applied')
    fl.image(s)
    #{{{alignment
    s.ift(['ph1','ph2'])
    s.ft('t2')
    opt_shift,sigma = correl_align(s,indirect_dim='nScans',
            ph1_selection = signal_pathway['ph1'],
            ph2_selection = signal_pathway['ph2'],sigma=50)
    s.ift('t2')
    s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename=None
    fl.next(r'after correlation, $\varphi$ domain')
    fl.image(s)
    s.ift('t2')
    s.ft(['ph1','ph2'])
    fl.next('after correl - time domain')
    fl.image(s)
    s.ft('t2')
    fl.next('after correl - freq domain')
    fl.image(s)
    s.ift('t2')
    s *= phasing
    s.ft('t2')
    fl.next('normal')
    fl.image(s)
    s = s['t2':(0,None)]
    s['t2':0] *= 0.5
    fl.next('processed data')
    fl.image(s)
    print(ndshape(s))
    s.reorder(['ph1','ph2','nScans','t2'])
    s.ift('t2')
    n_repeats = 32
    for j in range(n_repeats):
        data = s.C
        data /= sqrt(ndshape(data)['ph1']*ndshape(data)['ph2']) #normalization
        dt = diff(data.getaxis('t2')[r_[0,1]]).item()
        data.ft('t2')
        data /= sqrt(ndshape(data)['t2']) * dt
        manual_bounds = select_pathway(data['t2':(0,200)],signal_pathway)
        N = ndshape(manual_bounds)['t2']
        df = diff(data.getaxis('t2')[r_[0,1]]).item()
        manual_bounds.integrate('t2')
        all_results['repeats',j] = manual_bounds
    std_off_pathway = (
            data['ph1',0]['ph2',0]['t2':(0,200)].C.run(lambda x:
                abs(x)**2).mean_all_but(['t2']).mean('t2').run(sqrt))
    error_pathway = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    s_int,frq_slice = integral_w_errors(data,signal_pathway,error_pathway,
            indirect='nScans',fl=fl,return_frq_slice=True)
    fl.next('integrated for error')
    fl.plot(s_int,'.',capsize=6,alpha=0.3,label='integral_w_error')
    propagated_variance_from_inactive = N * df ** 2 *std_off_pathway **2
    manual_bounds.set_error(sqrt(propagated_variance_from_inactive.data.real))
    fl.plot(manual_bounds,'.',capsize=6,label=r'propagated from inactive std',alpha=0.5)
    all_results.mean('repeats',std=True)
    manual_bounds.set_error(all_results.get_error())
    fl.plot(manual_bounds,'.',capsize=6,label=r'manual std from repeats',alpha=0.5)
    #s_std = data.C.run(np.std,'t2')
    #manual_bounds.set_error(s_std)
    #fl.plot(manual_bounds,'.',capsize=6,label='numpy std',alpha=0.5)
    #fl.show();quit()
    fl.show();quit()



