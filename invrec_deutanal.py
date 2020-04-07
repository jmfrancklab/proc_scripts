print("Running")
from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from scipy.optimize import minimize,curve_fit,least_squares
from numpy import random
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
#baseline fitting
fl = figlist_var()
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
        ('w2_200309',2),
        ]:
    d = find_file(exp_name, exp_type='NMR_Data_AG', dimname = 'indirect', expno=expno)
print(ndshape(d))

print(d.getaxis('indirect'))
d.chunk('indirect',['indirect','ph1','ph2'],[-1,4,2]) #expands the indirect dimension into indirect, ph1, and ph2. inner most dimension is the inner most in the loop in pulse sequence, is the one on the farthest right. brackets with numbers are the number of phase cycle steps in each one. the number of steps is unknown in 'indirect' and is therefore -1.
print(d.getaxis('indirect'))
print(d.getaxis('ph1'))
print(d.getaxis('ph2'))
d.setaxis('ph1',r_[0:4.]/4) #setting values of axis ph1 to line up
d.setaxis('ph2',r_[0:2.]/4) #setting values of axis ph1 to line up
d.setaxis('indirect', d.get_prop('vd'))
fl.next('time domain') #switch to time domain as a string based name for fig
fl.image(d) #labeling
fl.next('FT + coherence domain')
#titling to coherence domain
d.ft('t2',shift=True) #fourier transform
d.ft(['ph1','ph2']) #fourier transforming from phase cycle dim to coherence dimension
d.reorder(['indirect','t2'], first=False)
print("after reorder",ndshape(d))
fl.image(d)
#fl.show();quit()
fl.next('select coherence pathway')
d = d['ph2',0]['ph1',0] # this brings you from 50 indirect dimensions to 5, The grouped rows represent ph1, we see signal in the top group which is 3 or -1, the signal is in the bottom of the two (which represent ph2) and so we select ph2, 0
fl.image(d)
d.reorder('t2')
print(ndshape(d))
#fl.show();quit()

# Use this to get the 1st order phase correction
SW = diff(d.getaxis('t2')[r_[0,-1]]).item()
thisph1 = nddata(r_[-6:6:2048j],'phi1')
oned_data = d['indirect',15].C
phase_test_r = oned_data*exp(-1j*2*pi*thisph1/SW*oned_data.fromaxis('t2'))
phase_test_rph0 = phase_test_r.C.sum('t2')
phase_test_rph0 /= abs(phase_test_rph0)
phase_test_r /= phase_test_rph0
cost = abs(phase_test_r.real).sum('t2')
fl.next('phasing cost function for first order correction')
fl.plot(cost,'.')
#fl.show()
#quit()
# Read the mininum from the above plot to get the first order correction
fl.next('phased first dimension')
fl.plot(oned_data,label='before')

# This phases the first indirect dimension (indirect 0)
ph1_0dim = -1.070
d['indirect',0] *= exp(-1j*2*pi*ph1_0dim/SW*d['indirect',0].fromaxis('t2')) # first order corr
# determining 0th order correction
d_ph0 = d['indirect',0].C.sum('t2') # summing along t2
d_ph0 /= abs(d_ph0) # dividing by abs val of the sum
d['indirect',0] /= d_ph0
fl.plot(d['indirect',0],label='after')
#read min from this then plug in for ph1_1dim and also go back to cost function section and change oned_data =d['indirect',0].C to oned_data =d['indirect',1]. run and plug in for ph1_1dim value 
#quit()
# This phases the second indirect dimensions (indirect 1)
ph1_1dim = -1.826
d['indirect',1] *= exp(-1j*2*pi*ph1_1dim/SW*d['indirect',1].fromaxis('t2')) # first order corr
# determining 0th order correction
d_ph0 = d['indirect',1].C.sum('t2') # summing along t2
d_ph0 /= abs(d_ph0) # dividing by abs val of the sum
d['indirect',1] /= d_ph0
#quit()

ph1_2dim = 2.924
d['indirect',2] *= exp(-1j*2*pi*ph1_2dim/SW*d['indirect',2].fromaxis('t2'))
d_ph0 = d['indirect',2].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',2] /= d_ph0
#quit()

ph1_3dim = -0.876
d['indirect',3] *= exp(-1j*2*pi*ph1_3dim/SW*d['indirect',3].fromaxis('t2'))
d_ph0 = d['indirect', 3].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',3] /= d_ph0
#quit()

ph1_4dim = 1.891
d['indirect',4] *= exp(-1j*2*pi*ph1_4dim/SW*d['indirect',4].fromaxis('t2'))
d_ph0 = d['indirect',4].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',4] /= d_ph0
#quit()

ph1_5dim = 5.729
d['indirect',5] *= exp(-1j*2*pi*ph1_5dim/SW*d['indirect',5].fromaxis('t2'))
d_ph0 = d['indirect',5].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',5] /= d_ph0
#quit()

ph1_6dim = -0.091
d['indirect',6] *= exp(-1j*2*pi*ph1_6dim/SW*d['indirect',6].fromaxis('t2'))
d_ph0 = d['indirect',6].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',6] /= d_ph0
#quit()

ph1_7dim = 1.959
d['indirect',7] *= exp(-1j*2*pi*ph1_7dim/SW*d['indirect',7].fromaxis('t2'))
d_ph0 = d['indirect',7].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',7] /= d_ph0
#quit()

ph1_8dim = -0.846
d['indirect',8] *= exp(-1j*2*pi*ph1_8dim/SW*d['indirect',8].fromaxis('t2'))
d_ph0 = d['indirect',8].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',8] /= d_ph0
#quit()

ph1_9dim = 4.670
d['indirect',9] *= exp(-1j*2*pi*ph1_9dim/SW*d['indirect',9].fromaxis('t2'))
d_ph0 = d['indirect',9].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',9] /= d_ph0
#quit()

ph1_10dim = -0.976
d['indirect',10] *= exp(-1j*2*pi*ph1_10dim/SW*d['indirect',10].fromaxis('t2'))
d_ph0 = d['indirect',10].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',10] /= d_ph0
#quit()

ph1_11dim = -0.812
d['indirect',11] *= exp(-1j*2*pi*ph1_11dim/SW*d['indirect',11].fromaxis('t2'))
d_ph0 = d['indirect',11].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',11] /= d_ph0
#quit()

ph1_12dim = -0.941
d['indirect',12] *= exp(-1j*2*pi*ph1_12dim/SW*d['indirect',12].fromaxis('t2'))
d_ph0 = d['indirect',12].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',12] /= d_ph0
#quit()

ph1_13dim = 5.695
d['indirect',13] *= exp(-1j*2*pi*ph1_13dim/SW*d['indirect',13].fromaxis('t2'))
d_ph0 = d['indirect',13].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',13] /= d_ph0
#quit()

ph1_14dim = -1.133
d['indirect',14] *= exp(-1j*2*pi*ph1_14dim/SW*d['indirect',14].fromaxis('t2'))
d_ph0 = d['indirect',14].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',14] /= d_ph0
#quit()

ph1_15dim = -1.046
d['indirect',15] *= exp(-1j*2*pi*ph1_15dim/SW*d['indirect',15].fromaxis('t2'))
d_ph0 = d['indirect',15].C.sum('t2')
d_ph0 /= abs(d_ph0)
d['indirect',15] /= d_ph0
#quit()
print(ndshape(d))
#for x in range(1,17):
 #   d['indirect',-1*x] *= -1.0
#d['indirect',0] *= -1
#d['indirect',1] *= -1
#d['indirect',2] *= -1
#d['indirect',3] *= -1
#d['indirect',4] *= -1
#d['indirect',5] *= -1
#d['indirect',6] *= -1
d['indirect',7] *= -1
d['indirect',8] *= -1
d['indirect',9] *= -1
d['indirect',10] *= -1
d['indirect',11] *= -1
d['indirect',12] *= -1
d['indirect',13] *= -1
d['indirect',14] *= -1
d['indirect',15] *= -1
fl.next('Plotting phased spectra')
for j in range(ndshape(d)['indirect']):
    fl.plot(d['indirect',j]['t2':(-50,50)],
        alpha=0.5,
        label='vd=%g'%d.getaxis('indirect')[j])
#exponential curve
rec_curve = d['t2':(-50,50)].C.sum('t2')
fl.next('recovery curve')
fl.plot(rec_curve,'o')
#fl.show();quit()

#estimating T1
min_index = abs(d).run(sum, 't2').argmin('indirect',raw_index=True).data
min_vd = d.getaxis('indirect')[min_index]
est_T1 = min_vd/log(2)
print("Estimated T1 is:", est_T1,"s")


#attempting ILT plot with NNLS_Tikhonov_190104

T1 = nddata(logspace(-5,1,150),'T1')
l = sqrt(logspace(-3.0,0.5,35)) 

def vec_lcurve(l):
    return d.C.nnls('indirect',T1,lambda x,y: 1.0-2*exp(-x/y), l=l)
x=vec_lcurve(l) 
x_norm = x.get_prop('nnls_residual').data
r_norm = x.C.run(linalg.norm,'T1').data

fl.next('L-Curve')
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
#fl.show();quit()
this_l = 0.045 #pick number in l curve right before it curves up
soln = d.real.C.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
soln.reorder('t2',first=False)
soln.rename('T1','log(T1)')
soln.setaxis('log(T1)',log10(T1.data))
fl.next('solution')
fl.image(soln)
fl.show();quit()
print("SAVING FILE")
np.savez('ag_'+exp_name+'_'+str(expno)+'_ILT_inv',
        data=soln.data,
        logT1=soln.getaxis('log(T1)'),
        t2=soln.getaxis('t2'))               
print("FILE SAVED")
quit()

