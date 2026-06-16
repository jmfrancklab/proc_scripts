"""Convert synthetic 2D data to tables of integrals.

This is the updated version of ``examples/broken/generate_integrals.py``.
The old helper named ``generate_integrals`` is no longer exported by
``pyspecProcScripts``; the current supported route is
``rough_table_of_integrals``.
"""

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp


np.random.seed(2021)
plt.rcParams["image.aspect"] = "auto"

t2, vd, power, ph1, ph2 = sp.symbols("t2 vd power ph1 ph2")
echo_time = 10e-3

examples = [
    dict(
        label="IR",
        expression=23
        * (1 - 2 * sp.exp(-vd / 0.2))
        * sp.exp(1j * 2 * sp.pi * 100 * t2 - abs(t2) * 50 * sp.pi),
        variable_defs=OrderedDict(
            [
                ("vd", psd.nddata(np.r_[0:1:40j], "vd")),
                ("ph1", psd.nddata(np.r_[0:4] / 4.0, "ph1")),
                ("ph2", psd.nddata(np.r_[0, 2] / 4.0, "ph2")),
                ("t2", psd.nddata(np.r_[0:0.2:256j] - echo_time, "t2")),
            ]
        ),
        signal_pathway={"ph1": 0, "ph2": 1},
        indirect="vd",
        signal_range=(-400, 400),
    ),
    dict(
        label="Enhancement",
        expression=23
        * (1 - (32 * power / (0.25 + power)) * 150e-6 * 659.33)
        * sp.exp(1j * 2 * sp.pi * 100 * t2 - abs(t2) * 50 * sp.pi),
        variable_defs=OrderedDict(
            [
                ("power", psd.nddata(np.r_[0:4:25j], "power")),
                ("ph1", psd.nddata(np.r_[0:4] / 4.0, "ph1")),
                ("t2", psd.nddata(np.r_[0:0.2:256j] - echo_time, "t2")),
            ]
        ),
        signal_pathway={"ph1": 1},
        indirect="power",
        signal_range=(-200, 600),
    ),
]


with psd.figlist_var() as fl:
    for cfg in examples:
        fl.basename = f"({cfg['label']})"
        data = psd.fake_data(
            cfg["expression"],
            cfg["variable_defs"],
            cfg["signal_pathway"],
        )
        data.reorder([cfg["indirect"], "t2"], first=False)
        data.ft("t2")
        data /= np.sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        data_int, ax = prscr.rough_table_of_integrals(
            data,
            signal_range=cfg["signal_range"],
            signal_pathway=cfg["signal_pathway"],
            fl=fl,
            title=cfg["label"],
        )
        fl.next(f"{cfg['label']} table of integrals")
        fl.plot(data_int, "o")
