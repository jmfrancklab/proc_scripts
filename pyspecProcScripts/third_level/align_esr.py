from ..first_level.QESR_rescale import QESR_scalefactor
from collections import OrderedDict
import numpy as np
from numpy import r_, pi
import matplotlib.pyplot as plt


def align_esr(data_dict, correlation_slice=None, on_crossing=False, fl=None):
    r"""open and correlation align ESR spectra.
    Store the shifts and scalefactors as
    properties of the nddata, which are returned
    as a dictionary of nddata.

    Details
    =======

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

    Parameters
    ==========
    correlation_slice: tuple default None
        This should ONLY be needed for VERY
        different spectra -- e.g. pure MTSL vs.
        a protein -- where without it, the
        different hyperfines show up out of
        register.
        You can look at the correlation plot,
        and select a range of shifts for it to
        find a max over.
    """
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
    for label_str, d in data_dict.items():
        # {{{ load, rescale
        d.setaxis(Bname, lambda x: x / 1e4).set_units(Bname, "T")
        d /= QESR_scalefactor(d)
        if "harmonic" in d.dimlabels:
            d = d["harmonic", 0]
        # set up so that we FT into a symmetric time domain
        d.set_ft_initial(Bname, "f")
        d -= d[Bname, -50:].C.mean(Bname).data
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
    BSW = maxB - minB
    # }}}
    # {{{ arrange the figures in the PDF
    if fl:
        fl.par_break()  # each fig on new line
        fl.next("Raw")
        fl.par_break()  # each fig on new line
        fl.next("correlation", legend=True)
        fl.par_break()
        fl.next(
            "aligned, autoscaled", figsize=(3 * 1.05 * 1.618, 3), legend=True
        )
        fl.par_break()
        fl.next("centered spectra", figsize=(3 * 1.05 * 1.618, 3), legend=True)
    # }}}
    # {{{ pull the reference (largest) up front
    all_files.move_to_end(ref_spec, last=False)
    # }}}
    for j, (label_str, d) in enumerate(
        (k, v) for (k, v) in all_files.items() if v != ref_spec
    ):
        # {{{ just show the raw data
        if fl:
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
        else:
            normfactor = d.data.max()
            d.ift(Bname)
            correlation = d * ref_spec
            correlation /= normfactor  # just for display purposes, since only
            #                           the argmax is used
            correlation.ft_new_startpoint(
                Bname, "f"
            )  # I want to calculate things in terms of an offset,
            #    which can be positive or negative, and need to shift in
            #    the next step
            correlation.ft(Bname, shift=True)
            if fl:
                fl.next("correlation")
                fl.plot(correlation, label=label_str)
            if correlation_slice is not None:
                correlation = correlation[Bname:correlation_slice]
            thisshift = correlation.real.argmax(Bname).item()
            d *= np.exp(-1j * 2 * pi * d.fromaxis(Bname) * thisshift)
            d.ft(Bname)
            scaling = (d * ref_spec_Bdom).sum(Bname) / (
                ref_spec_Bdom * ref_spec_Bdom
            ).sum(Bname)
            scaling = scaling.real.item()
        if fl:
            fl.next("aligned, autoscaled")
        aligned_autoscaled[label_str] = d / scaling
        aligned_autoscaled[label_str].set_prop("scaling", scaling)
        if fl:
            fl.plot(
                aligned_autoscaled[label_str],
                label=f"{label_str}\nscaling {scaling}",
                alpha=0.5,
            )
    # {{{ this loop is to move all into the u domain and then find the average
    #     "center field"
    sum_abs = 0
    for label_str, d in aligned_autoscaled.items():
        if fl:
            fl.next("u domain")
        d.ift(Bname)
        if fl:
            fl.plot(
                d,
                label=f"{label_str}\nscaling {d.get_prop('scaling')}",
                alpha=0.5,
            )
            fl.plot(d)
        d.ft(Bname)
        sum_abs += abs(d)
    # }}}
    # {{{ pick the dip
    peak = sum_abs.C.argmax(Bname).item()
    if on_crossing:
        idx, dips = sum_abs.contiguous(
            lambda x: x / x.max() < 0.5, return_idx=True
        )
        dips = dips[np.where(np.diff(idx, axis=1).ravel() > 3)[0], :]
        center_point = (
            sum_abs[
                Bname : tuple(
                    dips[np.argmin(abs(np.mean(dips, axis=1) - peak)), :]
                )
            ]
            .argmin(Bname)
            .item()
        )
        if fl:
            fl.next("find center")
        sum_abs[Bname] -= center_point
        if fl:
            fl.plot(sum_abs, color="k")
    else:
        center_point = peak
    # }}}
    for label_str, d in aligned_autoscaled.items():
        d[Bname] -= center_point
        d.set_prop("Bcenter", center_point)
        if fl:
            fl.next("find center")
        if fl:
            fl.plot(abs(d))
        d.ift(Bname)
        if fl:
            fl.next("before centering -- ift")
        if fl:
            fl.plot(d)
        if fl:
            fl.next("after centering -- ift")
        d.ft_new_startpoint(Bname, "f", -BSW / 2, nearest=True)
        d.ft(Bname)
        d.ift(Bname)
        if fl:
            fl.plot(d)
            fl.next("centered spectra")
        d.ft(Bname)
        d = d[Bname : (-BSW / 2, BSW / 2)]
        d.human_units()
        Bname_new = f"$(B_0{-center_point/d.div_units(Bname,'T'):+#0.5g})$"
        d.rename(Bname, Bname_new)
        if fl:
            fl.plot(
                d,
                label=f"{label_str}\n×{d.get_prop('scaling')}",
                alpha=0.5,
                human_units=False,
            )
        aligned_autoscaled[label_str] = d
    if fl:
        fl.next("centered spectra")
        fl.adjust_spines("bottom")
        plt.title("")
        plt.ylabel("")
        plt.gca().set_yticks([])
    return aligned_autoscaled
