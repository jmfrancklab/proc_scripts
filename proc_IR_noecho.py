from pyspecdata import *
import numpy as np
from numpy import * 
from proc_scripts import *
from proc_scripts.fitting import recovery
from proc_scripts.load_data import postproc_dict
from sympy import symbols
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
#baseline fitting
t2 = symbols('t2')
logger = init_logging("info")
fl = fl_mod()
#loading data in:
#    this_l is the regularization lambda specific to
#    each dataset, which must be chosen from the
#    L-curve.
for searchstr,exp_type,which_exp,postproc,this_l,f_range in [
        ('w8_200224','test_equip',2,'ab_ir2h',0.008,(-150,150)),
        #('w12_200224',2,(-150,150)),
        #('ag_oct182019_w0_10',3,(-150,150)),
        #('ag_oct182019_w0_8',3,(-150,150)),
        #('ag_oct182019_w0_6',2,(-150,150)),
        #('ag_oct182019_w0_3',2,(-150,150)),
        #('ag_oct92019_w0_12',2,(-150,150)),
        #('ag_oct92019_w0_10',2,(-150,150)),
        #('ag_oct92019_w0_8',2,(-150,150)),
        #('ag_bulk_d20',2,(-150,150)),
        #('ag_sep062019_w0_12_IR',2,(-150,150)),
        #('ag_sep232019_w0_12_prod',2,(-150,150)),
        #('ag_sep232019_w0_8_prod',2,(-150,150)),
        #('ag_sep232019_w0_6_prod',2,(-150,150)),
        #('ag_sep232019_w0_3_prod',2,(-150,150)),
        #('ag_sep232019_w0_1_prod',2,(-150,150)),
        ]:
    s = find_file(searchstr, exp_type=exp_type, 
            expno=which_exp, 
            postproc=postproc, lookup=postproc_dict, 
            dimname='indirect')
    s.ift('t2')
    #{{{ select appropriate coherence channel
    s = s['ph2',0]['ph1',-1]
    s.reorder('t2')
    ph0 = zeroth_order_ph(s['t2':0],fl=fl)
    ph0 /= abs(ph0)
    s /= ph0
    #}}}
    s.ft('t2')
    #{{{exponential curve with fit
    f,T1 = recovery(s, f_range)
    fl.plot_curve(f,"inversion recovery curve")
    #}}}
    #{{{estimating T1
    min_index = abs(s).run(sum, 't2').argmin('indirect',raw_index=True).data
    min_vd = s.getaxis('indirect')[min_index]
    est_T1 = min_vd/log(2)
    print("Estimated T1 is:", est_T1,"s")
    #}}}
    #{{{attempting ILT plot with NNLS_Tikhonov_190104
    T1 = nddata(logspace(-3,1,150),'T1')
    l = sqrt(logspace(-8.0,0.05,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
    if this_l is None:
        def vec_lcurve(l):
            return s.nnls('indirect',T1,lambda x,y: 1.0-2*exp(-x/y), l=l)

        x=vec_lcurve(l) 

        x_norm = x.get_prop('nnls_residual').data
        r_norm = x.run(linalg.norm,'T1').data

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
    #}}}
    #{{{setting axis to incorporate SFO1 based off of 
    # acqu file of data
    sfo1 = 272.05
    arbitrary_reference = s.get_prop('acq')['BF1'] 
    logger.info(strm("SFO1 is",sfo1))
    s.setaxis('t2',lambda x:x + sfo1 - arbitrary_reference)
    #}}}
    #{{{creating plot off of solution to L curve
    if this_l is not None:
        soln = s.real.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
        soln.reorder('t2',first=False)
        soln.rename('T1','log(T1)')
        soln.setaxis('log(T1)',log10(T1.data))
        fl.next('solution')
        fl.image(soln['t2':(100,300)])
fl.show()
    #}}}
