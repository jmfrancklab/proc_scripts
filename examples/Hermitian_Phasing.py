r"""
Phasing and Timing Correction With Fake Data
============================================

Take fake data with a relatively symmetric echo
(:math:`T_2^*=1/50\pi`, echo time of 10 ms),
and demonstrate how we can automatically find the zeroth order phase and the
center of the echo in order to get data that's purely real in the frequency
domain.
"""

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd
import pyspecProcScripts as pypcs
import sympy as s

psd.init_logging(level="debug")

np.random.seed(2021)
plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1
t2, vd, power = s.symbols("t2 vd power")
# The fine sampling leaves enough acquired spectral width for the modeled
# dispersive tails, so this example does not need an artificial early slice.
synthetic_time_axis = np.r_[0:0.2:4096j] - 10e-3
with psd.figlist_var() as fl:
    for expression, orderedDict, signal_pathway, indirect, label in [
        (
            (
                23
                * (1 - 2 * s.exp(-vd / 0.2))
                * s.exp(+1j * 2 * s.pi * 100 * t2 - abs(t2) * 50 * s.pi)
            ),
            [
                ("vd", psd.nddata(np.r_[0:1:40j], "vd")),
                ("ph1", psd.nddata(np.r_[0:4] / 4.0, "ph1")),
                ("ph2", psd.nddata(np.r_[0, 2] / 4.0, "ph2")),
                ("t2", psd.nddata(synthetic_time_axis, "t2")),
            ],
            {"ph1": 0, "ph2": 1},
            "vd",
            "IR",
        ),
        (
            (
                23
                * (1 - (32 * power / (0.25 + power)) * 150e-6 * 659.33)
                * s.exp(+1j * 2 * s.pi * 100 * t2 - abs(t2) * 50 * s.pi)
            ),
            [
                ("power", psd.nddata(np.r_[0:4:25j], "power")),
                ("ph1", psd.nddata(np.r_[0:4] / 4.0, "ph1")),
                ("t2", psd.nddata(synthetic_time_axis, "t2")),
            ],
            {"ph1": 1},
            "power",
            "enhancement",
        ),
    ]:
        fl.basename = "(%s)" % label
        fig, ax_list = plt.subplots(1, 3, figsize=(7, 7))
        fig.suptitle(fl.basename)
        fl.next("Data processing", fig=fig)
        data = psd.fake_data(
            expression,
            OrderedDict(orderedDict),
            signal_pathway,
        )
        data.reorder([indirect, "t2"], first=False)
        data.ft("t2")
        data /= np.sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        fl.image(data, ax=ax_list[0])
        ax_list[0].set_title("Raw Data")
        # Keep the full acquired bandwidth so the homogeneous-linewidth fit
        # can determine how far the dispersive tails extend before slicing.
        data = pypcs.fid_from_echo(
            data,
            signal_pathway,
            max_alignment_shift=0,
            fl=fl,
        )
        print(
            label,
            "echo center:",
            data.get_prop("echo_center"),
            "inhomogeneous bounds:",
            data.get_prop("inh_bounds"),
            "homogeneous linewidth:",
            data.get_prop("homogeneous_linewidth"),
            "fallback:",
            data.get_prop("homogeneous_linewidth_is_fallback"),
            "conditioning rate:",
            data.get_prop("alignment_conditioning_rate"),
            "processing bounds:",
            data.get_prop("processing_bounds"),
        )
        fl.image(data, ax=ax_list[1], human_units=False)
        ax_list[1].set_title("Phased and centered (ν)")
        data.ift("t2")
        fl.image(data, ax=ax_list[2], human_units=False)
        ax_list[2].set_title("Phased and centered (t)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
