from pylab import *
from pyspecdata import *
t2 = nddata(r_[0:1:1024j], 't2')
vd = nddata(r_[0:1:15j], 'vd')
# this generates fake data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 50 Hz
data = 21*(1-2*exp(-vd/0.2)) * exp(+1j*2*pi*100*t2-t2*100*pi)
data['t2':0] *= 0.5
data.add_noise(0.5)
data.ft('t2', shift=True)
plot(data, alpha=0.5)
figure()
# here I plot w/ manually chosen integration bounds:
plot(data['t2':(-400,400)].integrate('t2'), 'o')
# Now you need to run your code that automatically chooses integration bounds
# and also assigns error
# NEW CODE HERE
show()
