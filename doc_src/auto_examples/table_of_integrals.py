"""Convert 2D to Integrals with Errors
===================================

Take a 2D dataset and convert it to a table of integrals with errors, utilizing
all the bells and whistles (frequency and time selection, alignment, etc.)

Demonstrate on a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise, and with fluctuating field
(normally distributed field variation).
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import integrate_limits, integral_w_errors
from numpy.random import normal, seed
seed(2021)
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 4

# {{{ generate the fake data
init_logging(level="debug")
fl = figlist_var()
t2 = nddata(r_[0:0.2:256j], "t2")
vd = nddata(r_[0:1:40j], "vd")
ph1 = nddata(r_[0, 2] / 4.0, "ph1")
ph2 = nddata(r_[0:4] / 4.0, "ph2")
signal_pathway = {"ph1": 0, "ph2": 1}
# this generates fake clean_data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 10 Hz
echo_time = 5e-3
clean_data = 21*(1 - 2*exp(-vd / 0.2))*exp(+1j*2*pi*100*(t2-echo_time) - abs(t2-echo_time)*50*pi)
clean_data *= exp(signal_pathway["ph1"]*1j*2*pi*ph1)
clean_data *= exp(signal_pathway["ph2"]*1j*2*pi*ph2)
clean_data["t2":0] *= 0.5
# {{{ model frequency drift
sigma1 = 0.05 # distribution in frequency domain (smaller is "smoother")
sigma2 = 0.003 # distribution in frequency domain (smaller is "smoother")
scale = 100 # amplitude of frequency variation
frq_noise = normal(scale=scale,size=ndshape(clean_data)["ph1"] * ndshape(clean_data)["ph2"] * ndshape(clean_data)["vd"])
frq_noise = frq_noise + 1j*normal(scale=scale,size=frq_noise.size)
# {{{ control the spectral density of the shifts to be gaussian
frq_noise = nddata(frq_noise,[-1],['temp'])
N = ndshape(frq_noise)['temp']
frq_noise.setaxis('temp',-0.5+r_[0:N]/N).set_units('temp','cycperscan')
frq_noise_dens = 5*exp(-frq_noise.fromaxis('temp')**2/2/sigma2**2)
frq_noise_dens += exp(-frq_noise.fromaxis('temp')**2/2/sigma1**2)
frq_noise *= frq_noise_dens
fl.next('frq-noise density')
fl.plot(frq_noise)
frq_noise.ift('temp')
frq_noise /= sqrt(ndshape(frq_noise)['temp']) * frq_noise.get_ft_prop('temp','df') # normalization
fl.next('frq-noise time domain')
fl.plot(frq_noise)
# }}}
frq_noise = nddata(frq_noise.data.real,[-1,2,4],['vd','ph1','ph2'])
# }}}
fake_data_noise_std = 1.0
clean_data.reorder(["ph1", "ph2", "vd"])
data = clean_data.C
data.add_noise(fake_data_noise_std)
data *= exp(1j*2*pi*frq_noise*(data.fromaxis('t2')-echo_time))
# at this point, the fake data has been generated
data.ft(["ph1", "ph2"])
# {{{ usually, we don't use a unitary FT -- this makes it unitary
data /= 0.5 * 0.25  # the dt in the integral for both dims
data /= sqrt(ndshape(data)["ph1"] * ndshape(data)["ph2"])  # normalization
# }}}
dt = diff(data.getaxis("t2")[r_[0, 1]]).item()
fl.next("fake data -- time domain")
fl.image(data)
data.ft("t2", shift=True)
# {{{ make it unitary again
data /= sqrt(ndshape(data)["t2"]) * dt
# }}}
fl.next("fake data -- freq domain")
fl.image(data)
# }}}
fl.show()