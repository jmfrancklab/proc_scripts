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
clock_correction=False
# }}}
for thisfile,exp_type,nodename,postproc,f_range,t_range,ILT in [
        ('210311_TEMPOL_500uM_cap_probe_33dBm','inv_rec','signal','spincore_IR_v1',
            (-0.119e3,0.225e3),(0,76e-3),False),
        #('w3_201111','test_equip',2,'ab_ir2h',(-200,200),(0,60e-3),False)
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    # {{{ since we need to relabel vd frequently, make a method
    def as_scan_nbr(d):
        return d.C.setaxis('vd','#').set_units('vd','scan #')
    # }}}
    s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
    s = s['t2':f_range]
    s.ift('t2')
    if 'indirect' in s.dimlabels:
        s.rename('indirect','vd')
    fl.next('time domain cropped log')
    fl.image(as_scan_nbr(s))
    # no rough centering anymore -- if anything, we should preproc based on τ,
    # etc, but otherwise let the hermitian test handle it
    #{{{centering, hermitian function test and zeroth order phasing
    rough_center = abs(s).convolve('t2',10).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.next('rough centering')
    fl.image(as_scan_nbr(s))
    #}}}
    if clock_correction:
        #{{{ clock correction
        clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
        s.ft('t2')
        fl.next('before clock correction')
        fl.image(as_scan_nbr(s))
        s_clock=s['ph1',0]['ph2',1].sum('t2')
        s.ift(['ph1','ph2'])
        min_index = abs(s_clock).argmin('vd',raw_index=True).item()
        s_clock *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
        s_clock['vd',:min_index+1] *=-1
        s_clock.sum('vd').run(abs)
        fl.next('clock correction')
        fl.plot(s_clock,'.',alpha=0.7)
        clock_corr = s_clock.argmax('clock_corr').item()
        pyplot.axvline(x=clock_corr, alpha=0.5, color='r')
        s *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
        fl.next('after auto-clock correction')
        s.ft(['ph1','ph2'])
        fl.image(s.C.setaxis('vd','#'))
        s.ift('t2')
        #}}}
    #{{{phasing the aligned data
    best_shift = hermitian_function_test(s['ph2',1]['ph1',0].C.mean('vd'))
    logger.info(strm("best shift is", best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('time domain after hermitian test')
    fl.image(as_scan_nbr(s).cropped_log())
    fl.next('frequency domain after hermitian test')
    s.ft('t2')
    fl.image(as_scan_nbr(s).cropped_log())
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
    fl.next('phased data frequency domain')
    s.ft('t2')
    fl.image(s.C.convolve('t2',10).setaxis('vd','#').set_units('vd','scan #'))
    fl.next('phased data freq domain')
    fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
    #}}}
    s.ift('t2')
    #}}}
    if postproc == 'spincore_IR_v1':
        #{{{subplot for imaging before alignment
        s.ft('t2')
        s.reorder('ph1',first=True)
        s.reorder('ph2',first=True)
        s.reorder('t2',first=False)
        fl.next('freq domain - before corr, overview')
        fl.image(s.C.setaxis(
    'vd','#').set_units('vd','scan #'))
        fig,(ax1,ax2) = plt.subplots(1,2)
        fl.next('before alignment', fig=fig)
        fl.image(s['ph2',ph2_val]['ph1',ph1_val].C.setaxis(
    'vd','#').set_units('vd','scan #'),ax=ax1)
        fl.image(abs(s)['ph2',ph2_val]['ph1',ph1_val].C.setaxis(
    'vd','#').set_units('vd','scan #'),ax=ax2)
        fl.next('time domain before alignment')
        s.ift('t2')
        fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
        s.ft('t2')
        #}}}
        #{{{ pre-alignment
        s_final = s.C
        s_before = s['vd',:3].C
        s_after = s['vd',3:].C
        s_after.ift(['ph2','ph1'])
        s_after.ift('t2')
        s_after.ft('t2')
        s_before.ift(['ph2','ph1'])
        s_before.ift('t2')
        s_before.ft('t2')
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
        s_after.ift('t2')
        s_after *= np.exp(-1j*2*pi*opt_shift_after*s_after.fromaxis('t2'))
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
        s_final['vd',:3] = s_before.C
        s_final['vd',3:] = s_after.C
        s = s_final.C
        fig,(ax1,ax2) = plt.subplots(1,2)
        fl.next('after alignment', fig=fig)
        fl.image(s['ph2',ph2_val]['ph1',ph1_val].C.setaxis(
    'vd','#').set_units('vd','scan #'),ax=ax1)
        fl.image(abs(s)['ph2',ph2_val]['ph1',ph1_val].C.setaxis(
    'vd','#').set_units('vd','scan #'),ax=ax2)
        fl.next('coherence domain-after corr')
        fl.image(as_scan_nbr(s))
        s.ift('t2')
        fl.next('time domain-after corr')
        fl.image(as_scan_nbr(s))
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
    fl.show()
    if ILT:
        T1 = nddata(logspace(-3,3,150),'T1')
        l = sqrt(logspace(-8.0,0.5,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
        plot_Lcurve =False
        if plot_Lcurve:
            def vec_lcurve(l):
                return s.C.nnls('indirect',T1,lambda x,y: 1.0-2*exp(-x/y), l=l)

            x=vec_lcurve(l) 

            x_norm = x.get_prop('nnls_residual').data
            r_norm = x.C.run(linalg.norm,'T1').data

            with figlist_var() as fl:
                fl.next('L-Curve')
                figure(figsize=(15,10))
                fl.plot(log10(r_norm[:,0]),log10(x_norm[:,0]),'.')
                annotate_plot = True
                show_lambda = True
                if annotate_plot:
                    if show_lambda:
                        for j,this_l in enumerate(l):
                            annotate('%0.3f'%this_l, (log10(r_norm[j,0]),log10(x_norm[j,0])),
                                    ha='left',va='bottom',rotation=45)
                    else:
                        for j,this_l in enumerate(l):
                            annotate('%d'%j, (log10(r_norm[j,0]),log10(x_norm[j,0])),
                                    ha='left',va='bottom',rotation=45)
            d_2d = s*nddata(r_[1,1,1],r'\Omega')
        offset = s.get_prop('proc')['OFFSET']
        this_l = 0.042#pick number in l curve right before it curves up
        o1 = 297.01 #o1 for free D2O
        sfo1 = s.get_prop('acq')['BF1']
        s.setaxis('t2',lambda x:
                x-o1)
        s.setaxis('t2',lambda x:
                x/(-sfo1)).set_units('t2','ppm')
        soln = s.real.C.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
        soln.reorder('t2',first=False)
        soln.rename('T1','log(T1)')
        soln.setaxis('log(T1)',log10(T1.data))
        fl.next('w=6 in hexane')
        fl.image(soln.C.set_units('t2','ppm'))
        logger.info(strm("SAVING FILE"))
        np.savez(searchstr+'_'+str(nodename)+'_ILT_inv',
                data=soln.data,
                logT1=soln.getaxis('log(T1)'),
                t2=soln.getaxis('t2'))               
        logger.info(strm("FILE SAVED"))

