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
from pyspecProcScripts import *
from numpy.random import normal, seed
import sympy as s
from collections import OrderedDict
seed(2021)
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 4

# {{{ generate the fake data
init_logging(level="debug")
fl = figlist_var()
# this generates fake clean_data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 10 Hz
t2, td, vd, ph1, ph2 = s.symbols('t2 td vd ph1 ph2')
echo_time = 5e-3
data = fake_data(
    21*(1 - 2*s.exp(-vd / 0.2))*s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi),
    OrderedDict([
        ("vd" , nddata(r_[0:1:40j], "vd")),
        ("ph1" , nddata(r_[0, 2] / 4.0, "ph1")),
        ("ph2" , nddata(r_[0:4] / 4.0, "ph2")),
        ("t2" , nddata(r_[0:0.2:256j]-echo_time, "t2"))]),
        {"ph1": 0, "ph2": 1})
data.reorder(["ph1", "ph2", "vd"])
fl.next("fake data -- time domain")
fl.image(data)
data.ft("t2")
# {{{ make it unitary again
data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop('t2','dt')
# }}}
fl.next("fake data -- freq domain")
fl.image(data)
# }}}
# AG: implement your function on this data -- then also modify the symbolic expression to get enhancement, etc, and do that as well
fl.show()
