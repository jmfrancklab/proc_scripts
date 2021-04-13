from pylab import *
from pyspecdata import *
from proc_scripts import integrate_limits, integral_w_errors
fl=figlist_var()
t2 = nddata(r_[0:1:1024j], 't2')
vd = nddata(r_[0:1:15j], 'vd')
ph1 = nddata(r_[0,2]/4.,'ph1')
ph2 = nddata(r_[0:4]/4.,'ph2')
signal_pathway = {'ph1':0,'ph2':1}
excluded_pathways = [(0,0),(0,3)]
# this generates fake data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 50 Hz
data = 21*(1-2*exp(-vd/0.2)) * exp(+1j*2*pi*100*t2-t2*100*pi)
data *= exp(1j*2*pi*ph1)
data *= exp(1j*2*pi*ph2)
data.ft(['ph1','ph2'])
data['t2':0] *= 0.5
data.add_noise(0.5)
data.ft('t2', shift=True)
data.reorder(['ph2','ph1','vd'])
fl.next('before integration')
fl.image(data, alpha=0.5)
figure()
# here I plot w/ manually chosen integration bounds:
plot(data['t2':(-400,400)].integrate('t2'), 'o')
# Now you need to run your code that automatically chooses integration bounds
# and also assigns error
# NEW CODE HERE
error_pathway = (set(((j,k) for j in range(2) for k in range(4)))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
data = integral_w_errors(data,signal_pathway,error_pathway)
fl.next('with error')
fl.plot(data,'o',label='real')
fl.plot(data.imag,'o',label='imaginary')

show()
