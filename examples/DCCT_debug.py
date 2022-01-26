from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from numpy.random import normal, seed
import sympy as s
from collections import OrderedDict

seed(2021)
# {{{ generate the fake data
init_logging(level="debug")
fl = fl_mod()
t2, td, vd,ph1, ph2 = s.symbols("t2 td vd ph1 ph2")
signal_pathway = {"ph1": 0, "ph2": 1}
t_range = (0, 40e-3)
echo_time = 10e-3
with figlist_var() as fl:
    for (
        expression,
        orderedDict,
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
            True,
            "IR",
            (-400, 400),
        )    ]:
        fl.basename = "(%s)" % label
        data = fake_data(expression, OrderedDict(orderedDict), signal_pathway)
        data.reorder([indirect, "t2"], first=False)
        data.ft("t2")
        # {{{ make it unitary again
        data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        # }}}
        myslice = data["t2":f_range]
        mysgn = determine_sign(select_pathway(myslice, signal_pathway))
        DCCT(myslice,figure(figsize=(20,10)),total_spacing = 0.1,
                LHS_pad = 0.25,RHS_pad=0.51,
                allow_for_text_default = 5,
                allow_for_ticks_default = 50,
                text_height=50,
                label_factor_default = 13,
                plot_title = 'Raw Data time')


