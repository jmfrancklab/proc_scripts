r"""
EPR correlation alignment
=========================

While we can align by microwave frequency and normalize
according to peak-to-peak amplitude, it scan still be
hard to identify subtle differences between ESR
spectra, and small imperfections -- such as free MTSL
-- can play an outsized role.

Therefore, here, we use correlation to align the
spectra and then use "dot-product scaling" to normalize
them.  By "dot-product scaling" we mean scaling the
amplitude of one vector (here a spectrum,
:math:`\mathbf{b}`) relative to a reference (here
:math:`\mathbf{a}`) to minimize the residual between the
two -- *i.e.* we minimize the expression

.. math::

    |\mathbf{a}-c\mathbf{b}|^2

by varying the
scaling constant :math:`c`.
The solution to this is

.. math::

    c = \frac{\Re[\mathbf{a}\cdot \mathbf{b}]}{||\mathbf{b}||^2}

In order to do all this, we need a common *x*-axis that
we can use for correlation, etc.
Here, we look for the fields that are furthest left and furthest right,
and for the smallest spacing between field samples
-- we use these values to construct a (therefore all-inclusive) x axis.

Also, for the purposes of dot-product scaling, 
it is better to scale the less noisy spectrum
(:math:`\mathbf{b}` above) 
relative to the noisier spectrum
(:math:`\mathbf{a}` above)
-- *i.e.* above, we want :math:`\mathbf{b}` to be less noisy.
Here, we simply find the largest spectrum in the group
(assuming it is least noisy) and use it as :math:`\mathbf{b}`.
"""

import pyspecdata as psd
from pyspecProcScripts import QESR_scalefactor
from collections import OrderedDict
import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
from numpy import r_, pi
from scipy.special import jv

mpl.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})


def check_startpoint(d):
    print(
        "check ft props:",
        {k: v for k, v in d.other_info.items() if k.startswith("FT")},
    )


psd.init_logging(level="debug")
filenames_w_labels = [
    ("220307_S175_KCl.DSC", "220307_S175_KCl"),
    ("220729_prS175.DSC", "220729 prS175"),
    ("220307_S175_KI.DSC", "220307_S175_KI"),
    ("220307_prS175_KH2PO4.DSC", "220307_S175_KH2PO4"),
]

Bname = "$B_0$"
all_files = OrderedDict()
aligned_autoscaled = {}
# {{{ load all files first and do the following:
#       -   determine ref_axis which spans all the axes with the finest
#           resolution → ref_axis
#       -   all_files is an ordered dict of all files
#       -   ref_spec, which is the label/key of the largest spectrum
all_axes_extents = []
maxval = 0
for j, (filename, label_str) in enumerate(filenames_w_labels):
    # {{{ load, rescale
    d = psd.find_file(filename, exp_type="francklab_esr/Farhana")
    d.setaxis(Bname, lambda x: x / 1e4).set_units(Bname, "T")
    d /= QESR_scalefactor(d)
    if "harmonic" in d.dimlabels:
        d = d["harmonic", 0]
    # set up so that we FT into a symmetric time domain
    d.set_ft_initial(Bname, "f")
    d -= d[Bname, :50].C.mean(Bname).data
    all_axes_extents.append(
        tuple(d.getaxis(Bname)[r_[0, -1]])  # min, max
        + (np.diff(d.getaxis(Bname)[r_[0, 1]]).item(),)  # difference
    )
    all_files[label_str] = d
    temp = d.data.max() - d.data.min()
    if temp > maxval:
        maxval = temp
        ref_spec = label_str
    # }}}
minB, maxB, dB = zip(*all_axes_extents)
minB = min(minB)
maxB = max(maxB)
dB = min(dB)
ref_axis = r_[minB : maxB + dB : dB]
BSW = maxB-minB
# }}}

with psd.figlist_var(width=0.7, filename="ESR_align_example.pdf") as fl:
    # {{{ arrange the figures in the PDF
    fl.par_break()  # each fig on new line
    fl.next("Raw")
    fl.par_break()  # each fig on new line
    fl.next("correlation", legend=True)
    fl.par_break()
    fl.next("aligned, autoscaled", figsize=(3 * 1.05 * 1.618, 3), legend=True)
    fl.par_break()
    fl.next("centered spectra", figsize=(3 * 1.05 * 1.618, 3), legend=True)
    # }}}
    # {{{ pull the reference (largest) up front
    all_files.move_to_end(ref_spec, last=False)
    # }}}
    all_shifts = []  # to store all the phase differences?
    for j, (label_str, d) in enumerate(
        (k, v) for (k, v) in all_files.items() if v != ref_spec
    ):
        # {{{ just show the raw data
        fl.next("Raw")
        fl.plot(d, label=label_str, alpha=0.5)
        # }}}
        d = d.interp(Bname, ref_axis.copy(), kind="linear")
        d.set_ft_initial(Bname, "f")
        if j == 0:
            ref_spec_Bdom = d.C
            ref_spec = d.C
            ref_spec.ift(Bname)
            ref_spec.run(np.conj)
            scaling = 1
            check_startpoint(ref_spec)
        else:
            normfactor = d.data.max()
            check_startpoint(ref_spec)
            d.ift(Bname)
            check_startpoint(ref_spec)
            correlation = d * ref_spec
            correlation /= normfactor  # just for display purposes, since only
            #                           the argmax is used
            correlation.ft_new_startpoint(
                Bname, 'f'
            )  # I want to calculate things in terms of an offset,
            #    which can be positive or negative, and need to shift in
            #    the next step
            correlation.ft(Bname, shift=True)
            fl.next("correlation")
            fl.plot(correlation, label=label_str)
            thisshift = correlation.real.argmax(Bname).item()
            d *= np.exp(-1j * 2 * pi * d.fromaxis(Bname) * thisshift)
            check_startpoint(d)
            d.ft(Bname)
            scaling = (d * ref_spec_Bdom).sum(Bname) / (
                ref_spec_Bdom * ref_spec_Bdom
            ).sum(Bname)
            scaling = scaling.real.item()
        fl.next("aligned, autoscaled")
        aligned_autoscaled[label_str] = d / scaling
        aligned_autoscaled[label_str].set_prop("scaling", scaling)
        fl.plot(
            aligned_autoscaled[label_str],
            label=f"{label_str}\nscaling {scaling}",
            alpha=0.5,
        )
    # {{{ this loop is to move all into the u domain and then find the average
    #     "center field"
    for label_str, d in aligned_autoscaled.items():
        fl.next("u domain")
        d.ift(Bname)
        fl.plot(
            d, label=f"{label_str}\nscaling {d.get_prop('scaling')}", alpha=0.5
        )
        fl.plot(d)
        d.ft(Bname)
        shift = abs(d).argmax(Bname).item()
        all_shifts.append(shift)
    # }}}
    mean_shift = np.mean(all_shifts)
    print(f"mean shift {mean_shift:#0.4g}")
    for label_str, d in aligned_autoscaled.items():
        d[Bname] -= mean_shift
        d.set_prop('Bcenter',mean_shift)
        d.ift(Bname)
        fl.next("before centering -- ift")
        fl.plot(d)
        fl.next("after centering -- ift")
        d.ft_new_startpoint(Bname, 'f', -BSW/2,
                            nearest=True)
        d.ft(Bname)
        d.ift(Bname)
        fl.plot(d)
        fl.next("centered spectra")
        d.ft(Bname)
        d = d[Bname:(-BSW/2,BSW/2)]
        d.human_units()
        Bname_new = f"$(B_0{-mean_shift/d.div_units(Bname,'T'):+#0.5g})$"
        d.rename(Bname,Bname_new)
        fl.plot(d, human_units=False)
    fl.next("centered spectra")
    # mpl.pyplot.legend("", frameon=False)
    # fl.autolegend_list[fl.current] = False
    fl.adjust_spines("bottom")
    # plt.savefig("single_mutant_overlay.pdf")
    plt.title("")
    plt.ylabel("")
    plt.gca().set_yticks([])
