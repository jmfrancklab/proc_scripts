from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import hermitian_function_test, zeroth_order_ph, recovery, integrate_limits, correl_align, ph1_real_Abs, postproc_dict
from sympy import symbols, latex, Symbol
from matplotlib import *
from scipy.signal import tukey
import matplotlib.pyplot as plt
import numpy as np
from sympy import exp as s_exp
init_logging('debug')
fl = figlist_var()
t2 = symbols('t2')
# {{{ input parameters
ph1_val = 0
ph2_val = 1
coh_err = {'ph1':1,# coherence pathways to use for error -- note that this
        #             should ideally be pathways that do NOT include any known
        #             artifacts
        'ph2':r_[0,2,3]}
# }}}
for thisfile,exp_type,nodename,postproc,f_range,t_range in [
        ('210127_OHTEMPO10mM_cap_probe_IR_1','test_equip','signal','spincore_IR_v1',
            (-0.5e3,0.5e3),(0,67e-3))
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
    s = s['t2':f_range]
    s.ift('t2')
    fl.next('raw data')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    for j in range(ndshape(s)['vd']):
        if s['vd',j].C.sum('t2')<1:
            s['vd',j] *= -1
    s_final = s.C * 0
    s.ft('t2')
    fl.next('before splitting and correl align')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    #{{{centering, hermitian function test and zeroth order phasing
    rough_center = abs(s).convolve('t2',10).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.next('rough centering')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    #}}}
    #{{{phasing the aligned data
    s.ift('t2')
    best_shift = hermitian_function_test(s['ph2',1]['ph1',0].C.mean('vd'))
    logger.info(strm("best shift is", best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('time domain after hermitian test')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #').cropped_log())
    fl.next('frequency domain after hermitian test')
    s.ft('t2')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #').cropped_log())
    s.ift('t2')
    ph0 = s['t2':0]['ph2',ph2_val]['ph1',ph1_val]
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('phased data time domain')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    s = s['t2':t_range]
    fl.next('phased data freq domain')
    s.ft('t2')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    #}}}
    #{{{slicing FID not sure if needed
    #s = s['t2':(0,None)]
    #s['t2',0] *= 0.5
    #fl.next('phased and FID sliced')
    #fl.image(s.C.setaxis(
    #'vd','#').set_units('vd','scan #'))
    #fl.next('phased and FID sliced -- frequency domain')
    #s.ft('t2')
    #fl.image(s.C.setaxis(
    #'vd','#').set_units('vd','scan #'))
    #}}}

    s.ift('t2')
    fl.next('t domain for apod')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    #{{{time shift and apodization
    t0 = s['ph2',ph2_val]['ph1',ph1_val].C.mean_all_but('t2').argmax('t2').item()
    s.setaxis('t2', lambda x: x - t0)
    apod = s.fromaxis('t2') * 0
    apod['t2':(None,30e-3)] = apod['t2':(None,30e-3)].run(lambda x: tukey(len(x)))
    s *= apod
    #}}}
    #{{{subplot for imaging before alignment
    s.ft('t2')
    s.reorder('ph1',first=True)
    s.reorder('ph2',first=True)
    s.reorder('t2',first=False)
    fl.next('freq domain - before corr, overview')
    fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
    fig,(ax1,ax2) = plt.subplots(1,2)
    fl.next('filtered and apodized -- before alignment', fig=fig)
    fl.image(s['ph2',ph2_val]['ph1',ph1_val].C.setaxis(
'vd','#').set_units('vd','scan #'),ax=ax1)
    fl.image(abs(s)['ph2',ph2_val]['ph1',ph1_val].C.setaxis(
'vd','#').set_units('vd','scan #'),ax=ax2)
    #}}}
    #{{{ pre-alignment
    s_before = s['vd',:2]
    s_after = s['vd',3:]
    s_after.ift(['ph2','ph1'])
    s_after.ift('t2')
    s_after.ft('t2')
    fl.next('before alignment - FT')
    fl.image(s_after.C.setaxis(
'vd','#').set_units('vd','scan #'), interpolation='bilinear')
    s_before.ift(['ph2','ph1'])
    s_before.ift('t2')
    s_before.ft('t2')
    fl.next('prior to alignment')
    fl.image(s_before, interpolation = 'bilinear')
    #}}}
    #{{{alignment for before and after null
    fl.basename='first pass'
    opt_shift_before,sigma_before = correl_align(s_before,indirect_dim='vd',
            ph1_selection=ph1_val,ph2_selection=ph2_val,fl=fl)
    fl.basename = None
    fl.basename ='first pass after null'
    opt_shift_after,sigma_after = correl_align(s_after,indirect_dim='vd',
            ph1_selection=ph1_val,ph2_selection=ph2_val,fl=fl)
    fl.basename = None
    s_before.ift('t2')
    s_before *= np.exp(-1j*2*pi*opt_shift_before*s_before.fromaxis('t2'))
    fl.next('out of correl alignment before')
    fl.image(s_before)
    s_after.ift('t2')
    s_after *= np.exp(-1j*2*pi*opt_shift_after*s_after.fromaxis('t2'))
    fl.next('out of correl alignment after')
    fl.image(s_after.C.setaxis(
'vd','#').set_units('vd','scan #'))
    s_before.reorder('ph1',first=True)
    s_before.reorder('ph2',first=True)
    s_before.reorder('t2',first=False)
    s_before.ft(['ph1','ph2'])
    s_after.reorder('ph1',first=True)
    s_after.reorder('ph2',first=True)
    s_after.reorder('t2',first=False)
    s_after.ft(['ph1','ph2'])
    #}}}
    #{{{recombining before and after null after alignment
    s_final['vd',:2]['t2',:68] = s_before
    s_final['vd',3:]['t2',103:] = s_after
    s = s_final
    fl.next('coherence domain-after corr')
    s.ft('t2')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    s.ift('t2')
    fl.next('time domain-after corr')
    fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    #}}}
    fl.next('recovery curve')
    s_signal = s['ph2',ph2_val]['ph1',ph1_val]
    # {{{ here we use the inactive coherence pathways to determine the error
    #     associated with the data
    s_forerror = s['ph2',coh_err['ph2']]['ph1',coh_err['ph1']]
    # variance along t2 gives error for the mean, then average it across all other dimension, then sqrt for stdev
    logger.info(strm(ndshape(s_forerror)))
    frq_slice = integrate_limits(s_signal)
    s_signal = s_signal['t2':frq_slice]
    s_forerror = s_forerror['t2':frq_slice]
    s_forerror.run(lambda x:abs(x)**2).mean_all_but(['vd','t2']).integrate('t2').run(sqrt)
    s_signal.integrate('t2').set_error(s_forerror.data)
    # }}}
    logger.info(strm("here is what the error looks like",s_signal.get_error()))
    fl.plot(s_signal,'o',label='real')
    fl.plot(s_signal.imag,'o',label='imaginary')
    fl.next('Spectrum - freq domain')
    s = s['ph2',ph2_val]['ph1',ph1_val]
    fl.plot(s)
    x = s_signal.fromaxis('vd')
    f = fitdata(s_signal)
    error = fitdata(s_forerror)
    M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd", real=True)
    error.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
    f.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
    error.fit()
    f.fit()
    logger.info(strm("output:",f.output()))
    logger.info(strm("latex:",f.latex()))
    T1 = 1./f.output('R_1')
    fl.next('fit')
    fl.plot(s_signal,'o', label='actual data')
    fl.plot(s_signal.imag,'o',label='actual imaginary')
    fl.plot(f.eval(100),label='fit')
    fl.plot(error.eval(100),label='fit error')
    plt.text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes,size='large',
            horizontalalignment='center',color='k')
    ax = plt.gca()
    fl.show();quit()

