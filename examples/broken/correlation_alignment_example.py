"""
Align data with significant frequency drift
===========================================

Takes a 2D data set and applies proper phasing corrections followed by 
aligning the data through a correlation routine.
"""

# TODO ☐: you need to make sure this is fixed so that it can be moved out of
#         the "broken" directory.
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
def frq_mask(s):
    """Generates a mask that is nonzero along frequencies only over the
    bandwidth of the signal
        Parameteres
        ===========
        s: nddata
            data that the mask is applied to

    Returns
        =======
        s: nddata
            copy of data with the mask applied
    """
    # we want to leave the original s unchanged and return a copy
    for_mask = s.C
    # {{{ find center frequency
    nu_center = psdpr.select_pathway(s.C.mean("repeats"),signal_pathway).C.argmax("t2")
    # }}}
    # {{{ Make mask using the center frequency and sigma (whose estimate here
    #     is 20)
    frq_mask = np.exp(
        -((for_mask.fromaxis("t2") - nu_center) ** 2) / (2 * 20.0**2)
    )
    # }}}
    for_mask.ift(list(signal_pathway))
    return for_mask * frq_mask


def Delta_p_mask(s, signal_pathway):
    """ Filters out all but the signal pathway and the "ph1":0 or
    {'ph1':0,'ph2':0} pathways (depending on which experiment below is used).
    Note this serves as an example function and other filter functions could
    alternatively be used"""
    for ph_name, ph_val in signal_pathway.items():
        s.ft(["Delta%s" % ph_name.capitalize()])
        s = (
            s["Delta" + ph_name.capitalize(), ph_val]
            + s["Delta" + ph_name.capitalize(), 0]
        )
    return s
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
        )
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
        mysgn = (
            psdpr.select_pathway(data, signal_pathway)
            .C.real.sum("t2")
            .run(np.sign)
        )
        #    this is the sign of the signal -- note how on the next line,
        #    I pass sign-flipped data, so that we don't need to worry about
        #    messing with the original signal
        data.ift(list(signal_pathway.keys()))
        opt_shift, sigma = psdpr.correl_align(
            data * mysgn,
            repeat_dims=indirect,
            signal_pathway=signal_pathway,
            sigma=3000 / 2.355,
            max_shift=300,  # this makes the Gaussian mask 3
            #                 kHz (so much wider than the signal), and
            #                 max_shift needs to be set just wide enough to
            #                 accommodate the drift in signal
            frq_mask_fn=frq_mask,
            ph_mask_fn=Delta_p_mask,
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
