"""
EPR correlation alignment
=========================

while we can align by microwave frequency and normalize
according to peak-to-peak amplitude, it scan still be
hard to identify subtle differences between ESR
spectra, and small imperfections -- such as free MTSL
-- can play an outsized role.

Therefore, here, we use correlation to align the
spectra and then use dot-product scaling to normalize
them.  Note that dot-product scaling gives the optimal
solution that minimizes the following expression by varying the
scaling constant :math:`c`:
:math:`|\mathbf{a}-c\mathbf{b}|^2`

In order to do all this, we need a common *x*-axis that
we can use for correlation, etc.
Here, we use the first spectrum in order to determine
the *x*-axis.
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import QESR_scalefactor
import matplotlib as mpl
import pickle
init_logging(level='debug')
filenames_w_labels =  [
        ('220307_S175_KCl.DSC','220307_S175_KCl'),
        ('220729_prS175.DSC', '220729 prS175'),
        ('220307_S175_KI.DSC','220307_S175_KI'),                
        ('220307_prS175_KH2PO4.DSC','220307_S175_KH2PO4'), 
B_name = "$B_0$"

with figlist_var(width=0.7, filename="ESR_align_example.pdf") as fl:
    # {{{ arrange the figures in the PDF
    fl.text(r"\par")  # each fig on new line
    fl.next("Raw")
    fl.text(r"\par")  # each fig on new line
    fl.next("correlation", legend=True)
    # }}}
    for j, (filename, label_str) in enumerate(filenames_w_labels):
        # {{{ load, rescale
        d = find_file(filename, exp_type="francklab_esr/Farhana")
        d /= QESR_scalefactor(d)
        if "harmonic" in d.dimlabels:
            d = d["harmonic", 0]
        d -= d[B_name, :50].C.mean(B_name).data
        # }}}
        # {{{ just show the raw data
        fl.next("Raw")
        fl.plot(d, label=label_str, alpha=0.5)
        # }}}
        if j == 0:
            # if the reference spectrum, store its range
            B_slice = d.getaxis(B_name)[r_[0, -1]]
        B_start, B_stop = (
            3400,
            3650,
        )  # assume this is a good range to cover any amount of Hall sensor drift
        ref_axis = r_[B_start:B_stop:4096j]
        d = d.interp(B_name, ref_axis.copy(), kind="linear")
        d[B_name] -= B_start
        if j == 0:
            # u_slice = (-3.41333333,3.41) # hard-coded based on smallest res
            ref_spec_Bdom = d.C
            ref_spec = d.C
            ref_spec.ft(B_name, shift=True)
            ref_spec.run(conj)
            scaling = 1
        else:
            normfactor = d.data.max()
            d.ft(B_name, shift=True)
            correlation = d * ref_spec
            correlation /= normfactor
            correlation.ft_clear_startpoints(
                B_name, t="reset"
            )  # because I want to calculate things in terms of an offset, which can be positive or negative, and need to shift in the next step
            correlation.ift(B_name, shift=True)
            fl.next("correlation")
            fl.plot(correlation, label=label_str)
            thisshift = correlation.real.argmax(B_name).item()
            d *= exp(1j * 2 * pi * d.fromaxis(B_name) * thisshift)
            d.ift(B_name)
            d.set_units(B_name, "G")
            scaling = (d * ref_spec_Bdom).sum(B_name) / (
                ref_spec_Bdom * ref_spec_Bdom
            ).sum(B_name)
            scaling = scaling.real.item()
        d[B_name] += B_start
        d = d[B_name:B_slice]
        fl.text(r"\par")  # each fig on new line
        fl.next("aligned, autoscaled", legend=True)
        fl.plot(d / scaling, label=f"{label_str}\nscaling {scaling}", alpha=0.5)
        fl.text(r"\par")  # each fig on new line
