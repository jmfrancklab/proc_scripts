"""
Phasing and Timing Correction With Fake Data
============================================

Take fake data with a relatively symmetric echo 
(:math:r`T_2^*=1/50$\pi$`, echo time of 10 ms),
and demonstrate how we can automatically find the zeroth order phase and the
center of the echo in order to get data that's purely real in the frequency
domain.
"""
import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt
import sympy as sp
from collections import OrderedDict
from numpy import random 
from numpy import r_,sqrt

random.seed(2021)
plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1
t2, td, vd, power, ph1, ph2 = sp.symbols("t2 td vd power ph1 ph2")
echo_time = 10e-3
f_range = (-400, 400)
with psd.figlist_var() as fl:
    for expression, orderedDict, signal_pathway, indirect, label in [
        (
            (
                23
                * (1 - 2 * sp.exp(-vd / 0.2))
                * sp.exp(+1j * 2 * sp.pi * 100 * t2 - abs(t2) * 50 * sp.pi)
            ),
            [
                ("vd", psd.nddata(r_[0:1:40j], "vd")),
                ("ph1", psd.nddata(r_[0:4] / 4.0, "ph1")),
                ("ph2", psd.nddata(r_[0, 2] / 4.0, "ph2")),
                ("t2", psd.nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 0, "ph2": 1},
            "vd",
            "IR",
        ),
        (
            (
                23
                * (1 - (32 * power / (0.25 + power)) * 150e-6 * 659.33)
                * sp.exp(+1j * 2 * sp.pi * 100 * t2 - abs(t2) * 50 * sp.pi)
            ),
            [
                ("power", psd.nddata(r_[0:4:25j], "power")),
                ("ph1", psd.nddata(r_[0:4] / 4.0, "ph1")),
                ("t2", psd.nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 1},
            "power",
            "enhancement",
        ),
    ]:
        fl.basename = "(%s)" % label
        fig, ax_list = plt.subplots(1, 4, figsize=(7, 7))
        fig.suptitle(fl.basename)
        fl.next("Data processing", fig=fig)
        data = psd.fake_data(expression, OrderedDict(orderedDict), signal_pathway)
        data.reorder([indirect, "t2"], first=False)
        data.ft("t2")
        data /= sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        fl.image(data, ax=ax_list[0])
        ax_list[0].set_title("Raw Data")
        data = data["t2":f_range]
        data.ift("t2")
        data /= prscr.zeroth_order_ph(prscr.select_pathway(data, signal_pathway), fl=fl)
        fl.image(data, ax=ax_list[1], human_units=False)
        ax_list[1].set_title("Zeroth Order \n Phase Corrected")
        fl.basename = "(%s)" % label
        best_shift = prscr.hermitian_function_test(
            prscr.select_pathway(data.C.mean(indirect), signal_pathway), fl=fl
        )
        data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        data.ft("t2")
        fl.image(data, ax=ax_list[2])
        ax_list[2].set_title("Hermitian Test (ν)")
        data.ift("t2")
        fl.image(data, ax=ax_list[3], human_units=False)
        ax_list[3].set_title("Hermitian Test (t)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
