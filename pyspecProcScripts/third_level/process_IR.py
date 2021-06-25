"""Analyze Inversion Recovery Data
==================================
Analyzes data acquired from an Inversion Recovery experiment.
"""
from pyspecdata import *
from sympy import exp as s_exp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, latex, Symbol
from pyspecProcScripts import *
t2 = symbols('t2')

def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
def as_scan_nbr(d):
    '''since we need to relabel vd frequently- we make a method'''
    return d.C.setaxis('vd','#').set_units('vd','scan #')
signal_pathway = {'ph1':0,'ph2':1}
excluded_pathways = [(0,0)]
def process_IR(s, this_l = 0.032,
        l = sqrt(np.logspace(-8.0,0.5,35)),
        clock_correction = True,
        W=6.2,
        f_range = (None,None),
        t_range = (None,83e-3),
        IR = True,
        flip=False,
        sgn = None,
        ILT=False,
        fl=None):
    """Fully processes inversion recovery experiments to produce a fit inversion
    recovery curve to extrapolate the T1 time from.

    Parameters
    ==========
    s:      nddata
    this_l: int
            knee of the L curve that gives best balance between signal and
    l:      ndarray
            spacing for the L curve minimization
    clock_correction:   boolean
                        especially needed at times in spincore to correct for 
                        drift in frequency and time.
    W:      int
            repetition delay set for FIR experiments
    f_range:    tuple
                range in which the signal resides in the frequency domain.
    t_range:    tuple
                range in which the signal resides in the time domain.
    IR:     boolean
            True for a true inversion recovery fitting, false for FIR fitting.
    flip:   boolean 
            At higher powers it's been noted the signal can flip signs sometimes
            inverting the recovery curve-this corrects for it
    sign:   nddata  
            the nddata with the collected signs of the data along the signal that
            upon multiplication will convert all signal to the same sign. Useful
            for datasets where a null is present.
    ILT:    boolean
            If true, will produce an ILT image plot for visualizing the T1 as a function
            of offset.
    fl:     boolean
            option to show figures produced by processing script.
        """        
    s *= sgn
    s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
    s.ift('t2')
    s.reorder(['ph1','ph2','vd','t2'])
    #{{{ Applying DC offset
    s.ift(['ph1','ph2'])
    t_start = t_range[-1] / 4
    t_start *= 3
    rx_offset_corr = s['t2':(t_start,None)]
    rx_offset_corr = rx_offset_corr.data.mean()
    s -= rx_offset_corr
    s.ft('t2')
    s.ft(['ph1','ph2'])
    #}}}
    zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('vd',raw_index=True).item()
    if 'indirect' in s.dimlabels:
        s.rename('indirect','vd')
    # no rough centering anymore -- if anything, we should preproc based on τ,
    # etc, but otherwise let the hermitian test handle it
    #{{{ phasing the data
    s = s['t2':f_range]
    s.ift('t2')
    if clock_correction:
        #{{{ clock correction
        clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
        s.ft('t2')
        if fl is not None:
            fl.next('before clock correction')
            fl.image(as_scan_nbr(s))
        s_clock=s['ph1',1]['ph2',0].sum('t2')
        s.ift(['ph1','ph2'])
        min_index = abs(s_clock).argmin('vd',raw_index=True).item()
        s_clock *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
        s_clock['vd',:min_index+1] *=-1
        s_clock.sum('vd').run(abs)
        if fl is not None:
            fl.next('clock correction')
            fl.plot(s_clock,'.',alpha=0.7)
        clock_corr = s_clock.argmax('clock_corr').item()
        plt.axvline(x=clock_corr, alpha=0.5, color='r')
        s *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
        s.ft(['ph1','ph2'])
        if fl is not None:
            fl.next('after auto-clock correction')
            fl.image(s.C.setaxis('vd','#'))
        s.ift('t2')
    #{{{Applying phase corrections    
    best_shift,max_shift = hermitian_function_test(select_pathway(s.C.mean('vd'),signal_pathway))
    logger.info(strm("best shift is", best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    if fl is not None:
        fl.next('time domain after hermitian test')
        fl.image(as_scan_nbr(s))
    s.ft('t2')
    if fl is not None:
        fl.next('frequency domain after hermitian test')
        fl.image(as_scan_nbr(s))
        #}}}
    s.ift('t2')
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    s /= ph0
    s.ft(['ph1','ph2'])
    if fl is not None:
        fl.next('zeroth order corrected')
        fl.image(as_scan_nbr(s))
    s.ft('t2')
    if fl is not None:
        fl.next('phased data -- frequency domain')
        fl.image(as_scan_nbr(s))
    #}}}
    #}}}
    if 'ph2' in s.dimlabels:
        s.reorder(['ph1','ph2','vd','t2'])
    else:
        s.reorder(['ph1','vd','t2'])
    #{{{Correlation Alignment
    fl.basename='correlation subroutine:'
    #for the following, should be modified so we can pass a mask, rather than specifying ph1 and ph2, as here
    s,opt_shift,sigma = correl_align(s,indirect_dim='vd',
            signal_pathway=signal_pathway,
            sigma=10)
    fl.basename = None
    if fl is not None:
        fl.next(r'after correlation, $\varphi$ domain')
        fl.image(as_scan_nbr(s))   
    s.ift('t2')
    s.ft(['ph1','ph2'])
    if fl is not None:
        fl.next(r'after correlation')
        fl.image(as_scan_nbr(s)) 
    if 'ph2' in s.dimlabels:
        s.reorder(['ph1','ph2','vd','t2'])
    else:
        s.reorder(['ph1','vd','t2']) 
    #}}}
    #{{{FID slice
    s = s['t2':(0,t_range[-1])]
    s['t2',0] *= 0.5
    s.ft('t2')
    if fl is not None:
        fl.next('FID sliced -- frequency domain')
        fl.image(as_scan_nbr(s))
    #}}}    
    s *= sgn
    data = s.C
    #zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('vd',raw_index=True).item()
    #if flip:
    #    s['vd',:zero_crossing] *= -1
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
    # }}}
    #{{{Integrating with associated error from excluded pathways    
    s_int,frq_slice = integral_w_errors(s,signal_pathway,error_path,
            fl=fl,return_frq_slice=True)
    x1 = s_int.get_error()
    x1[:] /= sqrt(2)
    logger.info(strm("here is what the error looks like",s_int.get_error()))
    if fl is not None:
        fl.next('Integrated data - recovery curve')
        fl.plot(s_int,'o',capsize=6, label='real')
        fl.plot(s_int.imag,'o',capsize=6,label='imaginary')    
    #}}}
    #{{{Fitting Routine
    x = s_int.fromaxis('vd')
    f = fitdata(s_int)
    M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd")
    if IR:
        logging.info(strm("fitting using the regular IR equation"))
        f.functional_form = Mi - 2*Mi*s_exp(-vd*R1)
    else:
        logging.info(strm("fitting using the FIR equation"))
        f.functional_form = Mi*(1-(2-s_exp(-W*R1))*s_exp(-vd*R1))
    f.fit()
    logger.info(strm("output:",f.output()))
    logger.info(strm("latex:",f.latex()))
    T1 = 1./f.output('R_1')
    if fl is not None:
        fl.next('fit',legend=True)
        fl.plot(s_int,'o', capsize=6, label='actual data')
        fl.plot(s_int.imag,'o',capsize=6,label='actual imaginary')
        fl.plot(f.eval(100),label='fit')
        plt.text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes,size='medium',
                horizontalalignment='center',verticalalignment='center',color='k',
                position=(0.33,0.95),fontweight='bold')
        plt.legend(bbox_to_anchor=(1,1.01),loc='upper left')
    print("YOUR T1 IS:",T1)
    return T1
    #}}}
    
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

