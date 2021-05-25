import pylab as plb
from pyspecdata import *
from scipy.optimize import minimize, leastsq
from sympy import exp as s_exp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, latex, Symbol
from proc_scripts import *
t2 = symbols('t2')

def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
def as_scan_nbr(d):
    '''since we need to relabel vd frequently- we make a method'''
    return d.C.setaxis('vd','#').set_units('vd','scan #')

def process_IR(s, label='', fl=None,
        this_l = 0.032,
        l = sqrt(np.logspace(-8.0,0.5,35)),
        signal_pathway = {'ph1':0,'ph2':1},
        excluded_pathways = [(0,0),(0,3)],
        coh_err = {'ph1':1,'ph2':r_[0,2,3]},
        clock_correction = True,
        W=6.2,
        f_range = (None,None),
        t_range = (None,83e-3),
        IR = True,
        ILT=False,
        plot_all=True
        ):
    s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
    s = s['t2':f_range]
    s.ift('t2')
    rx_offset_corr = s['t2':(0.02,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    if 'indirect' in s.dimlabels:
        s.rename('indirect','vd')
    fl.next('time domain')
    fl.image(as_scan_nbr(s))
    # no rough centering anymore -- if anything, we should preproc based on τ,
    # etc, but otherwise let the hermitian test handle it
    #{{{ phasing the aligned data
    best_shift = hermitian_function_test(s['ph2',1]['ph1',0].C.mean('vd'))
    logger.info(strm("best shift is", best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    if plot_all:
        fl.next('time domain after hermitian test')
        fl.image(as_scan_nbr(s))
    s.ft('t2')
    if plot_all:
        fl.next('frequency domain after hermitian test')
        fl.image(as_scan_nbr(s))
    s.ift('t2')
    if clock_correction:
        #{{{ clock correction
        clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
        s.ft('t2')
        if plot_all:
            fl.next('before clock correction')
            fl.image(as_scan_nbr(s))
        s_clock=s['ph1',0]['ph2',1].sum('t2')
        s.ift(['ph1','ph2'])
        min_index = abs(s_clock).argmin('vd',raw_index=True).item()
        s_clock *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
        s_clock['vd',:min_index+1] *=-1
        s_clock.sum('vd').run(abs)
        if plot_all:
            fl.next('clock correction')
            fl.plot(s_clock,'.',alpha=0.7)
        clock_corr = s_clock.argmax('clock_corr').item()
        plt.axvline(x=clock_corr, alpha=0.5, color='r')
        s *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
        s.ft(['ph1','ph2'])
        if plot_all:
            fl.next('after auto-clock correction')
            fl.image(s.C.setaxis('vd','#'))
        s.ift('t2')
        #}}}
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',0]['ph2',1] = 1
    phasing.ift(['ph1','ph2'])
    
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    s /= ph0
    s.ft(['ph1','ph2'])
    if plot_all:
        fl.next('zeroth order corrected')
        fl.image(s)
    s.ft('t2')
    if plot_all:
        fl.next('phased data -- frequency domain')
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
    if plot_all:
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
    if plot_all:
        fl.next(r'after correlation, $\varphi$ domain')
        fl.image(as_scan_nbr(s))
    s.ft(['ph1','ph2'])
    if plot_all:
        fl.next(r'after correlation')
        fl.image(as_scan_nbr(s))
    if 'ph2' in s.dimlabels:
        s.reorder(['ph1','ph2','vd','t2'])
    else:
        s.reorder(['ph1','vd','t2'])
    if plot_all:
        fl.next('after correlation -- frequency domain')
        fl.image(as_scan_nbr(s))
    s.ift('t2')
    if plot_all:
        fl.next('after correlation -- time domain')
        fl.image(as_scan_nbr(s))
    s = s['t2':(0,t_range[-1])]
    s['t2',0] *= 0.5
    if plot_all:
        # visualize time domain after filtering and phasing
        fl.next('FID sliced -- time domain')
        fl.image(as_scan_nbr(s))
    s.ft('t2')
    if plot_all:
        fl.next('FID sliced -- frequency domain')
        fl.image(as_scan_nbr(s))
    s['vd':(None,1)] *= -1
    # }}}
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
    # }}}
    s_signal = integral_w_errors(s,signal_pathway,error_path)
    logger.info(strm("here is what the error looks like",s_signal.get_error()))
    if plot_all:
        fl.next('Integrated data - recovery curve')
        fl.plot(s_signal,'o',label='real')
        fl.plot(s_signal.imag,'o',label='imaginary')
    s = select_pathway(s,signal_pathway)
    fl.next('Spectrum - freq domain')
    fl.plot(s)
    x = s_signal.fromaxis('vd')
    f = fitdata(s_signal)
    M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd")
    if IR:
        f.functional_form = Mi - 2*Mi*s_exp(-vd*R1)
    else:
        f.functional_form = Mi*(1-(2-s_exp(-W*R1))*s_exp(-vd*R1))
    f.fit()
    logger.info(strm("output:",f.output()))
    logger.info(strm("latex:",f.latex()))
    T1 = 1./f.output('R_1')

    fl.next('fit',legend=True)
    fl.plot(s_signal,'o', label='actual data')
    fl.plot(s_signal.imag,'o',label='actual imaginary')
    fl.plot(f.eval(100),label='fit')
    plt.text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes,size='medium',
            horizontalalignment='center',verticalalignment='center',color='k',
            position=(0.33,0.95),fontweight='bold')
    plt.legend(bbox_to_anchor=(1,1.01),loc='upper left')
    logger.info(strm("YOUR T1 IS:",T1))
    fl.show()
    return T1
    
    if ILT:
        T1 = nddata(np.logspace(-3,3,150),'T1')
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
        fl.image(soln)
        logger.info(strm("SAVING FILE"))
        if save_npz:
            np.savez(thisfile+'_'+str(nodename)+'_ILT_inv',
                    data=soln.data,
                    logT1=soln.getaxis('log(T1)'),
                    t2=soln.getaxis('t2'))               
        logger.info(strm("FILE SAVED"))
        T1_values[i] = T1
