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

# sphinx_gallery_thumbnail_number = 8

# {{{ generate the fake data
init_logging(level="debug")
fl = fl_mod()
t2, td, vd, power, ph1, ph2 = s.symbols('t2 td vd power ph1 ph2')
signal_pathway = {'ph1':0,'ph2':1}
t_range= (0,40e-3)
echo_time = 5e-3
with figlist_var() as fl:
    for expression, orderedDict, signal_pathway, indirect, label,f_range in [
        #(
        #    (
        #        21
        #        *(1 - 2*s.exp(-vd / 0.2))
        #        *s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi)
        #    ),
        #    [
        #        ("vd" , nddata(r_[0:1:40j], "vd")),
        #        ("ph1" , nddata(r_[0, 2] / 4.0, "ph1")),
        #        ("ph2" , nddata(r_[0:4] / 4.0, "ph2")),
        #        ("t2" , nddata(r_[0:0.2:256j]-echo_time, "t2"))
        #    ],
        #    {"ph1": 0, "ph2": 1},
        #    "vd",
        #    "IR",
        #    (-350,350)
        #),
        (
            (
                21
                * (1-(32*power/(0.25+power))*150e-6*659.33)
                *s.exp(+1j*2*s.pi*100*t2-abs(t2)*50*s.pi)
            ),
            [
                ("power",nddata(r_[0:4:25j],"power")),
                ("ph1",nddata(r_[0:4]/4.0,"ph1")),
                ("t2",nddata(r_[0:0.2:256j]-echo_time,"t2")),
            ],
            {"ph1":1},
            "power",
            "Enhancement",
            (-300,500)
        ),
        ]: 
        fl.basename = "(%s)"%label
        data = fake_data(expression, OrderedDict(orderedDict), signal_pathway)
        data.reorder([indirect,'t2'], first = False)
        data.ft("t2")
        # {{{ make it unitary again
        data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop('t2','dt')
        # }}}
        myslice = data['t2':f_range]
        mysgn = determine_sign(select_pathway(myslice, signal_pathway))
        data_int, data = process_data(s=data,signal_pathway=signal_pathway,
                searchstr=label, f_range=f_range,t_range=t_range, 
                sgn=mysgn, indirect=indirect,fl=fl)

