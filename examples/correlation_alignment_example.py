"""Align data with significant frequency drift
==============================================

Takes a 2D data set and applies proper phasing corrections followed by 
aligning the data through a correlation routine.
"""
from pyspecdata import *
from pyspecProcScripts import *
from pylab import *
import sympy as s
from collections import OrderedDict
from numpy.random import normal,seed
seed(2021)
rcParams['image.aspect'] = 'auto' #needed for sphinx gallery

#{{{generate fake data
fl = figlist_var()
#this generates fake clean data with a T2 of 0.2s
# amplitude of 23 (arbitrary), offset of 300 Hz,
# FWHM 10 Hz
f_range = (-400,400)
t2, td, vd, ph1, ph2 = s.symbols('t2 td vd ph1 ph2')
echo_time = 1e-3
data = fake_data(
        23*(1 - 2*s.exp(-vd / 0.2)) * s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi),
        OrderedDict([
            ('vd', nddata(r_[0:1:40j],'vd')),
            ('ph1', nddata(r_[0:4]/4.0,'ph1')),
            ('ph2', nddata(r_[0,2]/4.0,'ph2')),
            ('t2', nddata(r_[0:0.2:256j]-echo_time,'t2'))]),
            {'ph1':0,'ph2':1})
data.reorder(["ph1", "ph2", "vd"])
fl.next("fake data -- time domain")
fl.image(data)
data.ft("t2")
# {{{ make it unitary again
data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop('t2','dt')
# }}}
fl.next("fake data -- freq domain")
fl.image(data)
#}}}
#{{{Changing sign of all signal to be the same so we don't have to
# deal with the null of the signal
print(type(data))
myslice = data['t2':f_range]
mysgn = myslice['ph1',0]['ph2',1].real.sum('t2').run(np.sign)
data *= mysgn
#}}}
data = data['t2':f_range]
data.ift('t2')
#{{{Applying the phase corrections
best_shift,max_shift = hermitian_function_test(data['ph1',0]['ph2',1].C.mean('vd'))
data.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
fl.next('After hermitian function test -- Time domain')
fl.image(data)
data.ft('t2')
fl.next('After hermitian function test -- Frequency domain')
fl.image(data)
data.ift('t2')
ph0 = data['ph1',0]['ph2',1]['t2':0]
ph0 /= abs(ph0)
data /= ph0
fl.next('After phasing corrections applied')
fl.image(data)
#}}}
#{{{Applying Correlation Routine to Align Data
data.ft('t2')
data,opt_shift,sigma = correl_align(data,indirect_dim='vd',
        signal_pathway = {'ph1':0,'ph2':1},
        sigma = 50)
fl.next(r'After Correlation, $\varphi$ domain')
fl.image(data)
data.ift('t2')
data.ft(['ph1','ph2'])
fl.next('Aligned Data -- Time Domain')
fl.image(data)
data.ft('t2')
fl.next('Aligned Data -- Frequency Domain')
fl.image(data)
#}}}
fl.show();quit()



