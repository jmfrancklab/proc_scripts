from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import hermitian_function_test, zeroth_order_ph, recovery, integrate_limits, correl_align, ph1_real_Abs, postproc_dict,DCCT
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
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
signal_pathway = {'ph1':0,'ph2':1}
coh_err = {'ph1':1,# coherence pathways to use for error -- note that this
        #             should ideally be pathways that do NOT include any known
        #             artifacts
        'ph2':r_[0,2,3]}
# }}}
clock_correction=True
for thisfile,exp_type,nodename,postproc,f_range,t_range,IR,ILT in [
        ('210316_TEMPOL_1mM_cap_probe_34dBm','inv_rec','signal','spincore_IR_v1',
            (-0.256e3,-0.017e3),(0,44e-3),False,False),
        #('w3_201111','test_equip',2,'ag_IR2H',(-600,600),(0,None),True)
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    #{{{ since we need to relabel vd frequently- we make a method
    def as_scan_nbr(d):
        return d.C.setaxis('vd','#').set_units('vd','scan #')
    #}}}
    s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
    s = s['t2':f_range]
    s.ift('t2')
    if 'indirect' in s.dimlabels:
        s.rename('indirect','vd')
    fl.next('time domain')
    fl.image(as_scan_nbr(s))
    # no rough centering anymore -- if anything, we should preproc based on τ,
    # etc, but otherwise let the hermitian test handle it
    #{{{phasing the aligned data
    if nodename == 'signal':
        best_shift = hermitian_function_test(s['ph2',1]['ph1',0].C.mean('vd'))
        logger.info(strm("best shift is", best_shift))
        s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
        fl.next('time domain after hermitian test')
        fl.image(as_scan_nbr(s))
        fl.next('frequency domain after hermitian test')
        s.ft('t2')
        fl.image(as_scan_nbr(s))
        s.ift('t2')
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
    ph0 = select_pathway(s['t2':0], signal_pathway)
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('phased data -- frequency domain')
    s.ft('t2')
    fl.image(as_scan_nbr(s))
    #}}}
    #}}}
    if 'ph2' in s.dimlabels:
        s.reorder(['ph1','ph2','vd','t2'])
    else:
        s.reorder(['ph1','vd','t2'])
    zero_crossing = abs(select_pathway(s,signal_pathway)).sum('t2').argmin('vd', raw_index=True).item()
    logger.info(strm("zero crossing at",zero_crossing))
    s.ift(['ph1','ph2'])
    # {{{ so that the filter is in range, do a rough alignment
    frq_max = abs(s).argmax('t2')
    s.ift('t2')
    s *= np.exp(-1j*2*pi*frq_max*s.fromaxis('t2'))
    s.ft('t2')
    # }}}
    fl.next(r'after rough align, $\varphi$ domain')
    fl.image(as_scan_nbr(s))
    fl.basename='correlation subroutine -- before zero crossing:'
    #for the following, should be modified so we can pass a mask, rather than specifying ph1 and ph2, as here
    logger.info(strm("ndshape",ndshape(s),"zero crossing at",zero_crossing))
    if zero_crossing > 1:
        opt_shift,sigma = correl_align(s['vd',:zero_crossing+1],indirect_dim='vd',
                ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
                sigma=50)
        s.ift('t2')
        s['vd',:zero_crossing+1] *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
        s.ft('t2')
    else:
        logger.warning("You have 1 point or less before your zero crossing!!!!")
    fl.basename='correlation subroutine -- after zero crossing:'
    opt_shift,sigma = correl_align(s['vd',zero_crossing+1:],indirect_dim='vd',
            ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
            sigma=50)
    s.ift('t2')
    s['vd',zero_crossing+1:] *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename = None
    fl.next(r'after correlation, $\varphi$ domain')
    fl.image(as_scan_nbr(s))
    s.ft(['ph1','ph2'])
    fl.next(r'after correlation, DCCT')
    fl.image(as_scan_nbr(s))
    if 'ph2' in s.dimlabels:
        s.reorder(['ph1','ph2','vd','t2'])
    else:
        s.reorder(['ph1','vd','t2'])
    fl.next('after correlation -- frequency domain')
    fl.image(as_scan_nbr(s))
    s.ift('t2')
    fl.next('after correlation -- time domain')
    fl.image(as_scan_nbr(s))
    s = s['t2':(0,t_range[-1])]
    s['t2',0] *= 0.5
    # visualize time domain after filtering and phasing
    fl.next('FID sliced -- time domain')
    fl.image(as_scan_nbr(s))
    s.ft('t2')
    fl.next('FID sliced -- frequency domain')
    fl.image(as_scan_nbr(s))
    fl.next('Integrated data - recovery curve')
    s_signal = select_pathway(s,signal_pathway)
    # {{{ here we use the inactive coherence pathways to determine the error
    #     associated with the data
    s_forerror = s['ph2',coh_err['ph2']]['ph1',coh_err['ph1']]
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
    #fl.show();quit()
    fl.next('Spectrum - freq domain')
    s = select_pathway(s,signal_pathway)
    fl.plot(s)
    x = s_signal.fromaxis('vd')
    f = fitdata(s_signal)
    error = fitdata(s_forerror)
    M0,Mi,R1,vd,W = symbols("M_0 M_inf R_1 vd W", real=True)
    if IR:
        error.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
        f.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
    else:
        error.functional_form = Mi*(1-(2-s_exp(-W*R1))*s_exp(-vd*R1))
        f.functional_form = Mi*(1-(2-s_exp(-W*R1))*s_exp(-vd*R1))
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
    logger.info(strm("YOUR T1 IS:",T1))
    fl.show()
    if ILT:
        T1 = nddata(np.logspace(-3,3,150),'T1')
        l = sqrt(np.logspace(-8.0,0.5,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
        plot_Lcurve = False
        if plot_Lcurve:
            def vec_lcurve(l):
                return s.C.nnls('vd',T1,lambda x,y: 1.0-2*np.exp(-x/y), l=l)

            x=vec_lcurve(l) 

            x_norm = x.get_prop('nnls_residual').data
            r_norm = x.C.run(np.linalg.norm,'T1').data

            with figlist_var() as fl:
                fl.next('L-Curve')
                plt.figure(figsize=(15,10))
                fl.plot(np.log10(r_norm[:,0]),np.log10(x_norm[:,0]),'.')
                annotate_plot = True
                show_lambda = True
                if annotate_plot:
                    if show_lambda:
                        for j,this_l in enumerate(l):
                            plt.annotate('%0.3f'%this_l, (np.log10(r_norm[j,0]),np.log10(x_norm[j,0])),
                                    ha='left',va='bottom',rotation=45)
                    else:
                        for j,this_l in enumerate(l):
                            plt.annotate('%d'%j, (np.log10(r_norm[j,0]),np.log10(x_norm[j,0])),
                                    ha='left',va='bottom',rotation=45)
            d_2d = s*nddata(r_[1,1,1],r'\Omega')
        offset = s.get_prop('proc')['OFFSET']
        this_l = 0.032#pick number in l curve right before it curves up
        #o1 = 297.01 #o1 for free D2O
        o1 = s.get_prop('acq')['O1']
        sfo1 = s.get_prop('acq')['BF1']
        s.setaxis('t2',lambda x:
                x+o1)
        s.setaxis('t2',lambda x:
                x/(sfo1)).set_units('t2','ppm')
        s.set_prop('x_inverted',True)
        soln = s.real.C.nnls('vd',T1, lambda x,y: 1.0-2.*np.exp(-x/y),l=this_l)
        soln.reorder('t2',first=False)
        soln.rename('T1','log(T1)')
        soln.setaxis('log(T1)',np.log10(T1.data))
        fl.next('w=3')
        fl.image(soln)#.C.setaxis('t2','#').set_units('t2','ppm'))
        logger.info(strm("SAVING FILE"))
        #np.savez(thisfile+'_'+str(nodename)+'_ILT_inv',
        #        data=soln.data,
        #        logT1=soln.getaxis('log(T1)'),
        #        t2=soln.getaxis('t2'))               
        logger.info(strm("FILE SAVED"))
        fl.show()

