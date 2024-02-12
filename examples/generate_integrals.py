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
from numpy.random import seed
import sympy as s
from collections import OrderedDict

#init_logging(level="debug")

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 3

fl = figlist_var()#fl_mod()
t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
echo_time = 10e-3
with figlist_var() as fl:
    for (
        expression,
        variable_defs,
        signal_pathway,
        indirect,
        clock_correction,
        label,
        f_range,
    ) in [
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
            False,
            "IR",
            (-400, 400),
        ),
        (
            (
                23
                * (1 - (32 * power / (0.25 + power)) * 150e-6 * 659.33)
                * s.exp(+1j * 2 * s.pi * 100 * t2 - abs(t2) * 50 * s.pi)
            ),
            [
                ("power", nddata(r_[0:4:25j], "power")),
                ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
                ("t2", nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 1},
            "power",
            False,
            "Enhancement",
            (-200, 600),
        ),
    ]:
        fl.basename = "(%s)" % label
        data = fake_data(expression, OrderedDict(variable_defs), signal_pathway)
        data.ft("t2")
        # {{{ make data unitary again
        data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        # }}}
        data_int, data = generate_integrals(
            s=data,
            signal_pathway=signal_pathway,
            searchstr=label,
            f_range=f_range,
            indirect=indirect,
            clock_correction=clock_correction,
            fl=fl,
        )
