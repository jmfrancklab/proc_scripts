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
    s.reorder('t2',first=False)
    s.chunk('nScans',['repeats','nScans'],[32,-1])
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    fl.next('raw data')
    fl.image(s['nScans',2])
    s.ift('t2')
    fl.next('raw data time domain')
    fl.image(s['nScans',2])
    t_range=(0,0.06)
    f_range = (-1e3,1e3)
    all_results = ndshape(s)
    all_results.pop("t2").pop("ph1").pop("ph2")
    all_results = all_results.alloc()
    all_results.reorder(['nScans','repeats'])
    n_nScans = 32
    s.ift(['ph1','ph2'])
    print("2")
    rx_offset_corr = s['t2':(0.045,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    s = s['t2':f_range]
    print("3")
    s.ft(['ph1','ph2'])
    s.ift('t2')
    print("4")
    best_shift = hermitian_function_test(select_pathway(s['nScans',2],signal_pathway))
    print("5")
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',1]['ph2',0] = 1
    phasing.ift(['ph1','ph2'])
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
    s /= ph0
    s.ft(['ph1','ph2'])
    fl.next('phase corrections applied')
    fl.image(s)
    #{{{rough alignment
    s.ift(['ph1','ph2'])
    s.ft('t2')
    frq_max = abs(s).argmax('t2')
    s.ift('t2')
    s *= np.exp(-1j*2*pi*frq_max*s.fromaxis('t2'))
    s.ft('t2')
    #}}}
    fl.basename = 'correlation subroutine'
    opt_shift,sigma = correl_align(s['nScans',15],indirect_dim='repeats',
            ph1_selection = signal_pathway['ph1'],
            ph2_selection = signal_pathway['ph2'])
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
    s = s['t2':(0,None)]
    s['t2':0] *= 0.5
    for j in range(n_nScans):
        data = s.C
        data /= sqrt(ndshape(data)['ph1']*ndshape(data)['ph2'])
        dt = diff(data.getaxis('t2')[r_[0,1]]).item()
        data.ft('t2')
        data /= sqrt(ndshape(data)['t2']) * dt
        manual_bounds = select_pathway(data['t2':(0,200)]['repeats',j],signal_pathway)
        print(ndshape(manual_bounds))
        N = ndshape(manual_bounds)['t2']
        df = diff(data.getaxis('t2')[r_[0,1]]).item()
        manual_bounds.integrate('t2')
        print(ndshape(manual_bounds))
        all_results['repeats',j] = manual_bounds
    std_off_pathway = (
            data['ph1',0]['ph2',0]['t2':(0,200)].C.run(lambda x:
                abs(x)**2).mean_all_but(['t2','repeats'])
            .mean('t2').run(sqrt))
    print("off-pathway std", std_off_pathway / sqrt(2))
    propagated_variance_from_inactive = N * df ** 2 *std_off_pathway **2
    propagated_variance = N *df**2*2**2 *2
    fl.next('different types of error')
    manual_bounds.set_error(sqrt(propagated_variance))
    fl.plot(manual_bounds,".",capsize=6,label=r'propagated from programmed variance',
            alpha=0.5)
    manual_bounds.set_error(sqrt(propagated_variance_from_inactive.data))
    fl.plot(manual_bounds,'.',capsize=6,label=r'propagated from inactive std',alpha=0.5)
    all_results.mean('repeats',std=True)
    manual_bounds.set_error(all_results.get_error())
    fl.plot(manual_bounds,'.',capsize=6,label=r'std from repeats',alpha=0.5)
    fl.show();quit()
    error_pathway = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    s_int,frq_slice = integral_w_errors(s,signal_pathway,error_pathway,
            indirect='repeats',fl=fl,return_frq_slice=True)
    s = s_int
    fl.next('integrated for error')
    fl.plot(s,'o',capsize=6,alpha=0.3)
    fl.show();quit()



