from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from sympy import symbols
import os

fl = figlist_var()
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
d.reorder('indirect',first=False)
d.reorder('echo',first=True)
d.reorder('echo',first=False)
d.reorder('t2',first=True)
d.reorder('t2',first=False)
figure();title('request 4 (request 2 in time domain)')
image(d['ph8',0]['ph4',1]['m',1]['n',0])

s = d['ph8',0]['ph4',1]['m',1]['n',0].C

grad_list = r_[ 0.02, 0.15714286, 0.29428571, 0.43142857, 0.56857143, 0.70571429, 0.84285714, 0.98 ]
s.setaxis('indirect',0.535*grad_list)

s = s['indirect',0].C

echo_center = abs(s['echo',0]).argmax('t2').data.item()
s.setaxis('t2', lambda x: x-echo_center)
figure();title('time domain before hermitian test')
image(s)
best_shift,res = hermitian_function_test(s,fl=fl)
figure();title('residual from hermitian test')
plot(res)
show()
print("Best shift is",best_shift)
d.setaxis('t2', lambda x: x-echo_center-best_shift)
figure();title('time domain after hermitian test')
image(d)
figure();title('request 5 (after subtracting best shift)')
image(d['ph8',0]['ph4',1]['m',1]['n',0])
quit()
s.register_axis({'t2':0})
s = s['t2':(0,None)]
s.ft('t2')
figure();title('request 3')
image(s)
show();fl.show();quit()

f_axis = s.fromaxis('t2')
def costfun(p):
    zeroorder_rad,firstorder = p
    phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(-1j*2*pi*zeroorder_rad)
    corr_test = phshift * s
    return (abs(corr_test.data.imag)**2)[:].sum()
iteration = 0
def print_fun(x, f, accepted):
    global iteration
    iteration += 1
    logger.info(strm(iteration, x, f, int(accepted)))
    return
sol = basinhopping(costfun, r_[0.,0.],
        minimizer_kwargs={"method":'L-BFGS-B'},
        callback=print_fun,
        stepsize=100.,
        niter=100,
        T=1000.
        )
zeroorder_rad, firstorder = sol.x
phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
phshift *= exp(-1j*2*pi*zeroorder_rad)
s *= phshift
print("RELATIVE PHASE SHIFT WAS %0.1f\\us and %0.1f$^\circ$", firstorder, angle(zeroorder_rad)/pi*180)
figure();title('after phased - real ft')
image(s.real)
figure();title('after phased - imag ft')
image(s.imag)
if s['echo',0].data[:].sum().real < 0:
    s *= -1
figure();title('plot')
plot(s)
print("Ok")
show();quit()

np.savez('proc_DOSY_CPMG_'+filename+'_expno'+str(expno),
        data = final_spec.data,
        indirect = final_spec.getaxis('indirect'),
        echo = final_spec.getaxis('echo'),
        t2 = final_spec.getaxis('t2'),
        )
