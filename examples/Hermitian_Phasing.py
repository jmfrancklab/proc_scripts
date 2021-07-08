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
            {"ph1":1},
            "power",
            "enhancement",
        ),
    ]:
        fl.basename = label
        data = fake_data(expression, OrderedDict(orderedDict), signal_pathway)
        data.reorder([indirect, "t2"], first=False)
        data.ft("t2")
        data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        fl.next("Data in Frequency Domain")
        fl.image(data)
        data = data["t2":f_range]
        data.ift('t2')
        rough_center = abs(select_pathway(data,signal_pathway)).C.convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
        logger.info(strm('Rough center is:',rough_center))
        data.setaxis('t2', lambda x: x - rough_center).register_axis({"t2": 0})
        fl.next('Rough centering')
        data.ft('t2')
        mysgn = select_pathway(data, signal_pathway).C.real.sum("t2").run(np.sign)
        data *= mysgn
        fl.image(data)
        data.ift('t2')
        ph0 = select_pathway(data, signal_pathway)["t2":0]
        ph0 /= abs(ph0)
        data /= ph0
        fl.next("Zeroth order phasing correction applied")
        fl.image(data)
        #{{{ Applying the phase corrections
        best_shift, max_shift = hermitian_function_test(
            select_pathway(data.C.mean(indirect), signal_pathway),fl=fl)
        data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        data.ft('t2')
        data *= mysgn
        fl.next('After Hermitian Test, Frequency Domain')
        fl.image(data)
        data.ift("t2")
        fl.next('After Hermitian Test, Time domain')
        fl.image(data)
