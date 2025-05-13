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


# {{{ Define the frequency mask function and the ph cyc mask
def frq_mask(s, sigma=150.0):
    """Note that we assume that our mask is a product of a
    frequency-domain and a coherence-domain function.
    This returns a copy multiplied by the square root of the
    frequency-domain part,
    leaving the original data untouched.

    Parameters
    ==========
    s : nddata
        Signal, given in the frequency domain and coherence transfer
        (*vs.* phase) domain.
        The property `coherence_pathway` must be set.
    """
    assert s.get_ft_prop("t2")
    assert s.get_ft_prop(list(s.get_prop("coherence_pathway").keys())[0])
    # {{{ find center frequency
    nu_center = (
        psdpr.select_pathway(s, s.get_prop("coherence_pathway"))
        .mean("repeats")
        .argmax("t2")
    )
    # }}}
    # {{{ Make mask using the center frequency and sigma.
    #     Standard gaussian is 2σ² in the denominator -- the extra 2 is
    #     for sqrt.
    frq_mask = np.exp(
        -((s.fromaxis("t2") - nu_center) ** 2) / (4 * sigma**2)
    )
    # }}}
    # note that when we multiply, we automatically generate a copy
    return s * frq_mask


def Delta_p_mask(s):
    """Filters out all but the signal pathway and the "ph1":0 or
    {'ph1':0,'ph2':0} pathways (depending on which experiment below is used).
    Note this serves as an example function and other filter functions could
    alternatively be used"""
    for j, (ph_name, ph_val) in enumerate(
        s.get_prop("coherence_pathway").items()
    ):
        # TODO ☐ (later): leave for JF after all other todo's resolved and
        #         example runs:
        #
        #         Rather than actually applying the slice, return a
        #         boolean numpy array that can be used to slice the
        #         smooshed Δp dimension.
        #         This allows us to not only calculate the correlation
        #         function as a sum across Δp, but also to calculate the
        #         norm before such a sum.
        if j == 0:
            signal_path = s["Delta" + ph_name.capitalize(), ph_val]
            zero_path = s["Delta" + ph_name.capitalize(), 0]
        else:
            signal_path = signal_path["Delta" + ph_name.capitalize(), ph_val]
            zero_path = zero_path["Delta" + ph_name.capitalize(), 0]
    return signal_path + zero_path


# }}}
seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 4

t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
echo_time = 10e-3
f_range = (-400, 400)

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
        ).set_prop("coherence_pathway", signal_pathway)
        data.reorder([indirect, "t2"], first=False)
        data.ft("t2")
        data /= np.sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        psd.DCCT(  # note that fl.DCCT doesn't allow us to title the
            #        individual figures
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
        data["t2"] -= data.getaxis("t2")[0]  # needed for Hermitian Function
        #                                     (fid_from_echo does this
        #                                     automatically)
        best_shift = psdpr.hermitian_function_test(
            psdpr.select_pathway(data.C.mean(indirect), signal_pathway)
        )
        data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        data.ft("t2")
        psd.DCCT(data, bbox=gs[1], fig=fig, title="Phased and \n Centered")
        # }}}
        # {{{ Applying Correlation Routine to Align Data
        mysgn = (  # this is the sign of the signal -- note how on the next
            #        line, I pass sign-flipped data, so that we don't need to
            #        worry about messing with the original signal
            psdpr.select_pathway(data, signal_pathway)
            .C.real.sum("t2")
            .run(np.sign)
        )
        opt_shift = psdpr.correl_align(
            data * mysgn,
            frq_mask_fn=frq_mask,
            Delta_p_mask_fn=Delta_p_mask,
            repeat_dims=indirect,
            signal_pathway=signal_pathway,
            max_shift=300,  # this makes the Gaussian mask 3
            #                 kHz (so much wider than the signal), and
            #                 max_shift needs to be set just wide enough to
            #                 accommodate the drift in signal
            fl=fl,
        )
        # The shift correction should be made in the time and phase cycling
        # domain
        # (this is how it is done in the function itself as well)
        data.ift("t2").ift(list(signal_pathway))
        data *= np.exp(-1j * 2 * np.pi * opt_shift * data.fromaxis("t2"))
        data.ft("t2").ft(
            list(signal_pathway.keys())
        )  # FT both domains for plotting
        psd.DCCT(data, bbox=gs[2], fig=fig, title=r"Aligned Data ($\nu$)")
        data.ift("t2")
        psd.DCCT(data, bbox=gs[3], fig=fig, title=r"Aligned Data ($t$)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        # }}}
