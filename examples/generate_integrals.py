"""Show synthetic integral/alignment examples.

This is the fit_T1 version of the fake-data section from the fit_ODNP branch.
It generates synthetic inversion-recovery and enhancement data, then compares
the DCCT-style ``table_of_integrals`` path against
``rough_table_of_integrals`` on the same fake data.
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
        center_aligned_peak=True,
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
        center_aligned_peak=False,
    ),
]


with psd.figlist_var() as fl:
    for cfg in examples:
        fl.basename = f"(fake {cfg['label']})"
        # {{{ Generate fake data and scale the direct FT like real data
        data = psd.fake_data(
            cfg["expression"],
            cfg["variable_defs"],
            cfg["signal_pathway"],
        )
        data.reorder([cfg["indirect"], "t2"], first=False)
        data.set_prop("acq_params", {"tau_us": echo_time * 1e6})
        data.ft("t2")
        data /= np.sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        # }}}

        # {{{ Compare correlation alignment and rough alignment integrals
        corr_int, _ = prscr.table_of_integrals(
            data.C,
            signal_range=cfg["signal_range"],
            signal_pathway=cfg["signal_pathway"],
            fl=fl,
            title=f"fake {cfg['label']} correlation alignment",
            repeat_dims=cfg["indirect"],
            center_aligned_peak=cfg["center_aligned_peak"],
            propagate_error=False,
        )
        rough_int, _ = prscr.rough_table_of_integrals(
            data.C,
            signal_range=cfg["signal_range"],
            signal_pathway=cfg["signal_pathway"],
            fl=fl,
            title=f"fake {cfg['label']} rough alignment",
        )
        # }}}

        # {{{ Overlay the two integral tables
        fl.next(f"fake {cfg['label']} table-of-integrals comparison")
        fl.plot(
            corr_int.C.set_plot_color("k"),
            "o",
            label="fake table_of_integrals correlation alignment",
        )
        fl.plot(
            rough_int.C.set_plot_color("r"),
            "x",
            label="fake rough_table_of_integrals alignment",
        )
        plt.gca().set_title(
            f"fake {cfg['label']} table-of-integrals comparison"
        )
        plt.legend()
        # }}}
