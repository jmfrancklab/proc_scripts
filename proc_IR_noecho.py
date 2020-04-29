print("Running")
from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from scipy.optimize import minimize,curve_fit,least_squares
from numpy import random
from proc_scripts import hermitian_function_test,zeroth_order_ph, load_data_bruker
from sympy import symbols
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
#baseline fitting
fl = figlist_var()
t2 = symbols('t2')
def calc_baseline(this_d,
        ph1lim,
        npts=5,
        guess=None,
        show_plots=True):
    if show_plots: fl.next('try baseline correction')
    if show_plots: fl.plot(this_d,
            label='before')
    this_d_tdom = this_d.C.ift('t2')
    blank_tdom = this_d_tdom.C
    blank_tdom.data[:] = 0
    def vec_to_params(ini_vec):
        phcorr0,phcorr1 = ini_vec[:2]
        baseline_vec = ini_vec[2:].view(complex128)
        return phcorr0,phcorr1,baseline_vec
    def apply_corr(ini_vec):
        phcorr0,phcorr1,baseline_vec = vec_to_params(ini_vec)
        d_test = this_d.C
        d_test *= exp(-1j*phcorr1*d.fromaxis('t2')-1j*phcorr0)
        d_test.ift('t2')
        retval = d_test['t2',0:len(baseline_vec)] + baseline_vec
        d_test['t2',0:len(baseline_vec)] = retval
        return d_test.ft('t2')
    def generate_baseline(baseline_vec):
        baseline_data = blank_tdom.C
        baseline_data['t2',0:len(baseline_vec)] += baseline_vec
        return baseline_data
    def costfun(ini_vec):
        d_test = apply_corr(ini_vec)
        return abs(d_test.real).sum('t2').data.item()
    max_val = abs(this_d_tdom.data).max()
    print(max_val)
    mybounds = r_[-max_val,max_val][newaxis,:]*ones(npts*2)[:,newaxis]
    print(shape(mybounds))
    print(mybounds)
    mybounds = r_[
            r_[-pi,pi,-ph1lim,ph1lim].reshape(-1,2),
            mybounds]
    if guess is None:
        guess = zeros(npts*2+2)
    else:
        guess = r_[guess[0].real,
                guess[1].real,
                guess[2:].view(float64)]
    res = minimize(costfun, (guess,),
            method='L-BFGS-B',
            bounds=mybounds,
            )
    phcorr0,phcorr1,baseline_vec = vec_to_params(res.x)
    baseline = generate_baseline(baseline_vec)
    if show_plots: fl.plot(this_d*exp(-1j*phcorr1*d.fromaxis('t2')-1j*phcorr0)+baseline.C.ft('t2'),
            label='after')
    return phcorr0,phcorr1,baseline
#loading data in
for exp_name,expno in [
        ('w8_200224',2),
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
    d = load_data_bruker(exp_name,'IR_noecho',expno) 
    #titling to coherence domain
    #rough_center = abs(d)['ph2',0]['ph1',0].convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    #d.setaxis(t2-rough_center)
    d.ft('t2',shift=True) #fourier transform

    fl.next('time domain (all $\\Delta p$)')
    d.ift('t2')
    fl.image(d)
    fl.next('frequency domain (all $\\Delta p$)')
    d.ft('t2',pad=4096)
    fl.image(d)
    #fl.show();quit()
    fl.next('select coherence pathway and convolve')
    d = d['ph2',0]['ph1',-1].C # this brings you from 50 indirect dimensions to 5, The grouped rows represent ph1, we see signal in the top group which is 3 or -1, the signal is in the bottom of the two (which represent ph2) and so we select ph2, 0
    d.convolve('t2',5)
    fl.image(d)
    #fl.show();quit()
    #d.ift('t2')
    d.reorder('t2')
    fl.next('after 0th order correction')
    ph0 = zeroth_order_ph(d['t2':0],fl=None)
    ph0 /= abs(ph0)
    d /= ph0
    fl.plot(d)
    d *= -1
    #fl.show();quit()

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
print("Estimated T1 is:", est_T1,"s")

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
print("SFO1 is",sfo1)
d.setaxis('t2',lambda x:x + sfo1 - arbitrary_reference)
this_l = 0.071#pick number in l curve right before it curves up
soln = d.real.C.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
soln.reorder('t2',first=False)
soln.rename('T1','log(T1)')
soln.setaxis('log(T1)',log10(T1.data))
fl.next('solution')
fl.image(soln['t2':(100,300)])


fl.show();quit()
print("SAVING FILE")
np.savez(exp_name+'_'+str(expno)+'_ILT_inv',
        data=soln.data,
        logT1=soln.getaxis('log(T1)'),
        t2=soln.getaxis('t2'))               
print("FILE SAVED")
quit()


