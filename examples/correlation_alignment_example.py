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
from numpy.random import normal, seed

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
echo_time = 10e-3
f_range = (-400, 400)

with figlist_var() as fl:
    for expression, orderedDict, signal_pathway, indirect, label in [
        (
            (
                23
                * (1 - 2 * s.exp(-vd / 0.2))
                * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi)
            ),
            [
                ("vd", nddata(r_[0:1:40j], "vd")),
                ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
                ("ph2", nddata(r_[0, 2] / 4.0, "ph2")),
                ("t2", nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 0, "ph2": 1},
            "vd",
            "IR",
        ),
        (
            (
                23
                * (1 - (32 * power / (0.25 + power)) * 150e-6 * 659.33)
                * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi)
            ),
            [
                ("power", nddata(r_[0:4:25j], "power")),
                ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
                ("t2", nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 1},
            "power",
            "enhancement",
        ),
    ]:
        fl.basename = label
        data = fake_data(expression, OrderedDict(orderedDict), signal_pathway)
        data.reorder([indirect, "t2"], first=False)
        fl.next("Raw Data in Time Domain")
        fl.image(data)
        data.ft("t2")
        data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        fl.next("Data in Frequency Domain")
        fl.image(data)
        myslice = data["t2":f_range]
        mysgn = select_pathway(myslice, signal_pathway).real.sum("t2").run(np.sign)
        data *= mysgn
        data = data["t2":f_range]
        data.ift("t2")
        #{{{ Applying the phase corrections
        best_shift, max_shift = hermitian_function_test(
            select_pathway(data, signal_pathway).C.mean(indirect)
        )
        data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        fl.next("After hermitian function test -- Time domain")
        fl.image(data)
        data.ft("t2")
        fl.next("After hermitian function test -- Frequency domain")
        # both of these plots labeled "After hermitian function test" now show
        # the sign of the signal bouncing back and forth -- why is that!!!????
        fl.image(data)
        data.ift("t2")
        ph0 = select_pathway(data, signal_pathway)["t2":0]
        ph0 /= abs(ph0)
        data /= ph0
        fl.next("After phasing corrections applied")
        fl.image(data)
        #}}}
        #{{{ Applying Correlation Routine to Align Data
        data.ft("t2")
        data, opt_shift, sigma = correl_align(
            data, indirect_dim=indirect, signal_pathway=signal_pathway, sigma=50
        )
        data.ift("t2")
        for k, v in signal_pathway.items():
            data.ft([k])
        fl.next("Aligned Data -- Time Domain")
        fl.image(data)
        data.ft("t2")
        fl.next("Aligned Data -- Frequency Domain")
        fl.image(data)
        #}}} 
