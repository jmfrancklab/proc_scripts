from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols, latex, Symbol
from matplotlib import *
from scipy.signal import tukey
import matplotlib.pyplot as plt
import numpy as np
from sympy import exp as s_exp
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
init_logging('debug')
fl = fl_mod()
t2 = symbols('t2')
# {{{ input parameters
this_l = 0.032 # pick number in l curve right before it curves up
l = sqrt(np.logspace(-8.0,0.5,35)) # play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
signal_pathway = {'ph1':0,'ph2':1}
excluded_pathways = [(0,0),(0,3)] # exclude ph1 ph2, since it usually has a receiver glitch, 0,-1 is FID-like
coh_err = {'ph1':1,# coherence pathways to use for error -- note that this
        #             should ideally be pathways that do NOT include any known
        #             artifacts
        'ph2':r_[0,2,3]}
save_npz = False
# }}}
clock_correction=True
W = 6#repetition delay used for FIR

for thisfile,exp_type,nodename,postproc,f_range,t_range,IR,ILT in [
        ('210512_water_cap_probe_IR_33dBm','inv_rec','signal','spincore_IR_v1',
            (-0.473e3,0.8e3),(None,83e-3),True,False),
        #('w3_201111','test_equip',2,'ag_IR2H',(-600,600),(0,None),True)
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    print(ndshape(s))
    #s.set_units('t2','kHz')
    #fl.side_by_side('freq',s,f_range)
    #fl.show();quit()
    #{{{ since we need to relabel vd frequently- we make a method
    def as_scan_nbr(d):
        return d.C.setaxis('vd','#').set_units('vd','scan #')
    #}}}
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
    #fl.show();quit()
    # no rough centering anymore -- if anything, we should preproc based on τ,
    # etc, but otherwise let the hermitian test handle it
    #{{{ phasing the aligned data
    print("EDIT")
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
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',0]['ph2',1] = 1
    phasing.ift(['ph1','ph2'])
    fl.next('zeroth order corrected')
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    s /= ph0
    s.ft(['ph1','ph2'])
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
    #fl.show();quit()
    s = s['t2':(0,t_range[-1])]
    s['t2',0] *= 0.5
    # visualize time domain after filtering and phasing
    fl.next('FID sliced -- time domain')
    fl.image(as_scan_nbr(s))
    s.ft('t2')
    fl.next('FID sliced -- frequency domain')
    fl.image(as_scan_nbr(s))
    #s *= -1  
    s['vd':(None,1)] *= -1
    fl.next('Integrated data - recovery curve')
    # }}}
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
    # }}}
    s_signal = integral_w_errors(s,signal_pathway,error_path)
    logger.info(strm("here is what the error looks like",s_signal.get_error()))
    fl.plot(s_signal,'o',label='real')
    fl.plot(s_signal.imag,'o',label='imaginary')
    fl.next('Spectrum - freq domain')
    s = select_pathway(s,signal_pathway)
    fl.plot(s)
    x = s_signal.fromaxis('vd')
    print(ndshape(s))
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
    fl.next('fit')
    fl.plot(s_signal,'o', label='actual data')
    fl.plot(s_signal.imag,'o',label='actual imaginary')
    fl.plot(f.eval(100),label='fit')
    plt.text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes,size='large',
            horizontalalignment='center',color='k')
    ax = plt.gca()
    logger.info(strm("YOUR T1 IS:",T1))
    print("YOUR T1 IS:",T1)
    fl.show()
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
        fl.show()

