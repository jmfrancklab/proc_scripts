import pylab as plb
from pyspecdata import *
from scipy.optimize import minimize, leastsq
from sympy import exp as s_exp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, latex, Symbol
from pyspecProcScripts import *
from pyspecProcScripts.third_level.process_data import proc_data
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
        sign = None,
        ILT=False,
        fl=None):
    s_int = proc_data(s,label='Inversion Recovery',indirect='vd',fl=fl,
            f_range=f_range,t_range=t_range,clock_correction=clock_correction,flip=flip,
            sign=sign)
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
        f.functional_form = Mi - 2*Mi*s_exp(-vd*R1)
    else:
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

