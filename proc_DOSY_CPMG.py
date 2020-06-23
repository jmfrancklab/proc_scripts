from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
import os

t2 = symbols('t2')

expno = 3
filename = 'ab_jun172020_w0_5'
d = find_file(filename,
              expno=expno,
             exp_type='NMR_Data_AAB')

figure();title('raw data')
image(d)

l22 = d.get_prop('acq')['L'][22]
l25 = d.get_prop('acq')['L'][25]
d12 = d.get_prop('acq')['D'][12]
d11 = d.get_prop('acq')['D'][11]
p1 = d.get_prop('acq')['P'][1]
dwdel1 = 3.5e-6
anavpt = 1024
dwdel2 = (anavpt*0.05e-6)/2
TD = 65536
quadrature_points = TD/2
num_points_per_echo = quadrature_points/l25
acq_time = dwdel2*num_points_per_echo*2
tau_extra = 20e-6
tau_pad = tau_extra-6e-6
tau_pad_start = tau_extra-dwdel1-6e-6
tau_pad_end = tau_extra-6e-6
tE = dwdel1 + 5e-6 + tau_pad_start + 1e-6 + num_points_per_echo*(dwdel2*2) + tau_pad_end

d.chunk('t2',['echo','t2'],[int(l25),-1])
d.chunk('indirect',['indirect','phcyc'],[int(l22),-1])
d.chunk('phcyc',['ph8','ph4','m','n'],[2,2,2,2])
d.setaxis('ph8',r_[0.,2.]/4)
d.setaxis('ph4',r_[0.,2.]/4)
d.setaxis('m',r_[0,2.]/4)
d.setaxis('n',r_[0,2.]/4)
d.ft(['ph8','ph4','m','n'])
d.ft('t2',shift=True)
figure();title('request 2')
image(d)

s = d['ph8',0]['ph4',1]['m',1]['n',0].C

figure();title('signal')
image(abs(s.C.convolve('t2',100.)))

grad_list = r_[ 0.02, 0.15714286, 0.29428571, 0.43142857, 0.56857143, 0.70571429, 0.84285714, 0.98 ]
s.setaxis('indirect',0.535*grad_list)

#{{{ phasing code
def calc_baseline(this_d,
        ph1lim,
        npts=5,
        guess=None,
        show_plots=True):
    if show_plots: figure();title('try baseline correction')
    if show_plots: plot(this_d,
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
        d_test *= exp(-1j*phcorr1*this_d.fromaxis('t2')-1j*phcorr0)
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
    if show_plots:
        plot(this_d*exp(-1j*phcorr1*this_d.fromaxis('t2')-1j*phcorr0)+baseline.C.ft('t2'),
            label='after')
        legend()
    return phcorr0,phcorr1,baseline
#}}}

figure();title('Before phasing or baseline, full')
plot(s.real['indirect',0]['echo',0])
plot(s.imag['indirect',0]['echo',0])


SW = diff(s.getaxis('t2')[r_[0,-1]]).item()
thisph1 = nddata(r_[-2:2:2000j],'phi1')

oned_data = s['indirect',0]['echo',0]
print("ndshape s",ndshape(oned_data))

oned_data.ift('t2')
rough_center = abs(oned_data).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
oned_data.setaxis(t2-rough_center)
figure();title('rough centered data')
plot(oned_data)
residual,best_shift = hermitian_function_test(oned_data)
figure();title('hermitian test')
plot(residual)
print("best shift is",best_shift)
oned_data.ft('t2')
oned_data *= exp(-1j*2*pi*best_shift*s.fromaxis('t2'))
oned_data.ift('t2')
figure();title('time domain after hermitian test')
plot(oned_data)
ph0 = oned_data['t2':0]
oned_data /= ph0
oned_data = oned_data['t2':(0,None)]
oned_data.ft('t2')
figure();title('time domain after hermitian test and phasing')
plot(oned_data)


show();quit()



print("Ok")
quit()

phase_test_array = oned_data*exp(
        -1j*2*pi*thisph1/SW*oned_data.fromaxis('t2'))
phase_test_array_ph0 = phase_test_array.C.sum('t2')
phase_test_array_ph0 /= abs(phase_test_array_ph0)
phase_test_array /= phase_test_array_ph0
cost = abs(phase_test_array.real).sum('t2')
figure();title('cost function')
plot(cost,'.')
ph1_guess = cost.argmin('phi1').data.item()*2*pi/SW

#{{{ show_phasing func
def show_phasing(this_d):
    plot(this_d.angle/pi,'.')
    abs_spec = abs(this_d)
    abs_spec /= abs_spec.data.max()
    plot(abs_spec-1.0,'r',alpha=0.5)
#}}}
guess_included = s['indirect',0]['echo',0] * exp(
    -1j*s['indirect',0]['echo',0].fromaxis('t2')*ph1_guess)
ph0_guess = angle((abs(guess_included.C.sum('t2'))
        /guess_included.C.sum('t2')).data.item())
guess_included *= exp(-1j*ph0_guess)

figure();title('show phasing -- initial guess')
show_phasing(guess_included)

s *= exp(-1j*ph1_guess*s.fromaxis('t2')
        -1j*ph0_guess)

phcorr0,phcorr1,baseline = calc_baseline(s['indirect',0]['echo',0],
        abs(ph1_guess),
        show_plots=True,
        )# use the guessed phase correction as the

figure();title('full spectrum - after phasing')
plot(s.real['indirect',0]['echo',0])
plot(s.imag['indirect',0]['echo',0])
figure();title('sliced spectrum - after phasing')
plot(s.real['indirect',0]['echo',0]['t2':(-1000,1200)],label='real')
plot(s.imag['indirect',0]['echo',0]['t2':(-1000,1200)],label='imag')
legend()

SW = diff(s.getaxis('t2')[r_[0,-1]]).item()
thisph1 = nddata(r_[-6:6:2000j],'phi1')

# units of cycles per SW
oned_data = s['indirect',0]['echo',0]
print("ndshape s",ndshape(s))

phase_test_array = oned_data*exp(
        -1j*2*pi*thisph1/SW*oned_data.fromaxis('t2'))
phase_test_array_ph0 = phase_test_array.C.sum('t2')
phase_test_array_ph0 /= abs(phase_test_array_ph0)
phase_test_array /= phase_test_array_ph0
cost = abs(phase_test_array.real).sum('t2')
figure();title('cost function - second pass')
plot(cost,'.')

ph1_guess = cost.argmin('phi1').data.item()*2*pi/SW
guess_included = s['indirect',0]['echo',0] * exp(
    -1j*s['indirect',0]['echo',0].fromaxis('t2')*ph1_guess)
ph0_guess = angle((abs(guess_included.C.sum('t2'))
        /guess_included.C.sum('t2')).data.item())
guess_included *= exp(-1j*ph0_guess)

figure();title('show phasing -- initial guess second pass')
show_phasing(guess_included)
s *= exp(-1j*ph1_guess*s.fromaxis('t2')
        -1j*ph0_guess)
phcorr0,phcorr1,baseline = calc_baseline(s['indirect',0]['echo',0],
        abs(ph1_guess),
        show_plots=False,
        )# use the guessed phase correction as the
        #                limit for the size of the further phase
        #                correction
final_spec = s.C
final_spec['indirect',0] *= exp(-1j*phcorr1*s.fromaxis('t2')-1j*phcorr0)
final_spec['indirect',0] += baseline.C.ft('t2')

figure();title('baseline and phase corr second pass')
plot(s['indirect',0]['echo',0],label='before')
plot(-final_spec['indirect',0]['echo',0],label='after')
legend()

scaling = abs(s['t2',0:10]).sum('t2')
scaling /= scaling['indirect',0]
startguess = r_[phcorr0,phcorr1,baseline.data[:5].copy()]

for j in range(1,ndshape(s)['indirect']):
    figure();title('baseline and phase corr for %d'%(j+1))
    plot(s['indirect',j]['echo',0],label='before')
    thisguess = startguess.copy()
    thisguess[2:] *= scaling['indirect',j]['echo',0].data.item()
    phcorr0,phcorr1,baseline = calc_baseline(s['indirect',j]['echo',0],
                                            abs(startguess[1]*2),
                                            guess=thisguess,
                                            show_plots=False)
    final_spec['indirect',j] *= exp(-1j*phcorr1*s.fromaxis('t2')-1j*phcorr0)
    final_spec['indirect',j] += baseline.C.ft('t2')
    plot(final_spec['indirect',j],label='after')

print(ndshape(final_spec))

for x in range(int(l22)):
    if final_spec['indirect',x].sum('t2') < 0:
        final_spec['indirect',x].data *= -1

np.savez('proc_DOSY_CPMG_'+filename+'_expno'+str(expno),
        data = final_spec.data,
        indirect = final_spec.getaxis('indirect'),
        echo = final_spec.getaxis('echo'),
        t2 = final_spec.getaxis('t2'),
        )

show();quit()
