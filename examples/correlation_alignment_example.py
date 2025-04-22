"""
Align data with significant frequency drift
===========================================

Takes a 2D data set and applies proper phasing corrections followed by 
aligning the data through a correlation routine.
"""

import pyspecdata as psd
from pyspecdata import r_
import numpy as np
import pyspecProcScripts as psdpr
from pylab import rcParams
import matplotlib.pyplot as plt
import sympy as s
from collections import OrderedDict
from numpy.random import seed

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 4

t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
echo_time = 10e-3
f_range = (-400, 400)
def frq_mask(s,signal_pathway, direct="t2", indirect = "repeats", sigma = 20.0):
    assert s.get_ft_prop(direct), "You must be in the frequency domain!" 
    for phnames in signal_pathway.keys():
        assert not s.get_ft_prop( phnames), (
            str(phnames) + "must NOT be in coherence domain!"
        )
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    s.ft(list(signal_pathway))
    # {{{ find center frequency
    for j in range(len(signal_keys)):
        signal = s[signal_keys[j], signal_values[j]].C
    nu_center = signal.mean(indirect).C.argmax(direct)
    # }}}
    # {{{ center and mask using sigma
    frq_mask = np.exp(
            -((s.fromaxis(direct)-nu_center)**2) / (2* sigma**2)
            )
    # }}} 
    s.ift(list(signal_pathway))
    masked_s = s*frq_mask
    return masked_s, frq_mask 

with psd.figlist_var() as fl:
    for expression, orderedDict, signal_pathway, indirect, label in [
        (
            (
                23
                * (1 - 2 * s.exp(-vd / 0.2))
                * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi)
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
                * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi)
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
        # {{{ equivalent of subplot
        fig = plt.figure(figsize=(11, 7))
        gs = plt.GridSpec(1, 4, figure=fig, wspace=0.4)
        # }}}
        fig.suptitle(fl.basename)
        fl.next("Data Processing", fig=fig)
        data = psd.fake_data(
            expression, OrderedDict(orderedDict), signal_pathway
        )
        data.set_prop("coherence_pathway",signal_pathway)
        data.reorder([indirect, "t2"], first=False)
        data.ft("t2")
        data /= np.sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        psd.DCCT(  # note that fl.DCCT doesn't allow us to title the individual
            #        figures
            data,
            bbox=gs[0],
            fig=fig,
            title="Raw Data",
        )
        data = data["t2":f_range]
        data.ift("t2")
        data /= psdpr.zeroth_order_ph(
            psdpr.select_pathway(data, signal_pathway)
        )
        # }}}
        # {{{ Applying the phase corrections
        data["t2"] -= data.getaxis("t2")[0] # needed for Hermitian Function
        best_shift = psdpr.hermitian_function_test(
            psdpr.select_pathway(data.C.mean(indirect), signal_pathway)
        )
        data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        data.ft("t2")
        psd.DCCT(data, bbox=gs[1], fig=fig, title="Phased and \n Centered")
        # }}}
        # {{{ Applying Correlation Routine to Align Data
        mysgn = (
            psdpr.select_pathway(data, signal_pathway)
            .C.real.sum("t2")
            .run(np.sign)
        )
        #    this is the sign of the signal -- note how on the next line,
        #    I pass sign-flipped data, so that we don't need to worry about
        #    messing with the original signal
        data.ift(list(signal_pathway.keys()))
        opt_shift = psdpr.correl_align(
            data * mysgn,
            repeat_dims=[indirect],
            signal_pathway=signal_pathway,
            max_shift=300,  # this makes the Gaussian mask 3
            #                 kHz (so much wider than the signal), and
            #                 max_shift needs to be set just wide enough to
            #                 accommodate the drift in signal
            frq_mask_fn = frq_mask,
            fl=fl,
        )
        # removed display of the mask (I think that's what it was)
        data.ift("t2")
        data *= np.exp(-1j * 2 * np.pi * opt_shift * data.fromaxis("t2"))
        data.ft(list(signal_pathway.keys()))
        data.ft("t2")
        psd.DCCT(data, bbox=gs[2], fig=fig, title=r"Aligned Data ($\nu$)")
        data.ift("t2")
        psd.DCCT(data, bbox=gs[3], fig=fig, title=r"Aligned Data ($t$)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        # }}}
