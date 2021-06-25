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
            ('ph1', nddata(r_[0,2]/4.0,'ph1')),
            ('ph2', nddata(r_[0:4]/4.0,'ph2')),
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
#{{{Changing sign of all signal to be the same so we don't have to
# deal with the null of the signal
myslice = s['t2':f_range]
mysgn = determine_sign(s['ph1',0]['ph2',1])
s *= mysgn
#}}}
s = s['t2':f_range]
s.ift('t2')
#{{{Applying the phase corrections
best_shift,max_shift = hermitian_function_test(s['ph1',0]['ph2',1].C.mean('vd'))
s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
fl.next('After hermitian function test -- Time domain')
fl.image(s)
s.ft('t2')
fl.next('After hermitian function test -- Frequency domain')
fl.image(s)
s.ift('t2')

fl.show();quit()



