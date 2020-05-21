print("Running")
from pyspecdata import *
from numpy import random
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from sympy import symbols
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
#baseline fitting
fl = figlist_var()
t2 = symbols('t2')
#loading data in
for searchstr,exp_type,which_exp,postproc in [
        ('w8_200224','test_equip',2,'ag_IR2H'),
        #('w12_200224',2),
        #('ag_oct182019_w0_10',3),
        #('ag_oct182019_w0_8',3),
        #('ag_oct182019_w0_6',2),
        #('ag_oct182019_w0_3',2),
        #('ag_oct92019_w0_12',2),
        #('ag_oct92019_w0_10',2),
        #('ag_oct92019_w0_8',2),
        #('ag_bulk_d20',2),
        #('ag_sep062019_w0_12_IR',2),
        #('ag_sep232019_w0_12_prod',2),
        #('ag_sep232019_w0_8_prod',2),
        #('ag_sep232019_w0_6_prod',2),
        #('ag_sep232019_w0_3_prod',2),
        #('ag_sep232019_w0_1_prod',2),
        ]:
    d = find_file(searchstr, exp_type=exp_type, 
            expno=which_exp, 
            postproc=None, lookup=postproc_dict, 
            dimname='indirect') 
    fl.next('time domain (all $\\Delta p$)')
    d.ift('t2')
    fl.image(d)
    fl.next('frequency domain (all $\\Delta p$)')
    d.ft('t2',pad=4096)
    fl.image(d)
    fl.next('select coherence pathway and convolve')
    d.convolve('t2',5)
    fl.image(d)
    d.reorder('t2')
    ph0 = zeroth_order_ph(d['t2':0],fl=None)
    ph0 /= abs(ph0)
    d /= ph0
    d *= -1

fl.next('Plotting phased spectra')
for j in range(ndshape(d)['indirect']):
    fl.plot(d['indirect',j]['t2':(-150,150)],
        alpha=0.5,
        label='vd=%g'%d.getaxis('indirect')[j])

#exponential curve
rec_curve = d['t2':(-150,150)].C.sum('t2')
fl.next('recovery curve')
fl.plot(rec_curve,'o')
#fl.show()
#quit()

#estimating T1
min_index = abs(d).run(sum, 't2').argmin('indirect',raw_index=True).data
min_vd = d.getaxis('indirect')[min_index]
est_T1 = min_vd/log(2)
logger.info(strm("Estimated T1 is:", est_T1,"s"))

#attempting ILT plot with NNLS_Tikhonov_190104

T1 = nddata(logspace(-3,1,150),'T1')
l = sqrt(logspace(-3.0,0.005,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
plot_Lcurve = False
if plot_Lcurve:
    def vec_lcurve(l):
        return d.C.nnls('indirect',T1,lambda x,y: 1.0-2*exp(-x/y), l=l)

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
    d_2d = d*nddata(r_[1,1,1],r'\Omega')
#fl.show()
#quit()
sfo1 = 273.76
arbitrary_reference = d.get_prop('acq')['BF1'] # will eventually be 
logger.info(strm("SFO1 is",sfo1))
d.setaxis('t2',lambda x:x + sfo1 - arbitrary_reference)
this_l = 0.178#pick number in l curve right before it curves up
soln = d.real.C.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
soln.reorder('t2',first=False)
soln.rename('T1','log(T1)')
soln.setaxis('log(T1)',log10(T1.data))
fl.next('solution')
fl.image(soln['t2':(100,300)])
fl.show();quit()



