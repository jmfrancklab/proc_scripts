r"""
EPR correlation alignment
=========================

while we can align by microwave frequency and normalize
according to peak-to-peak amplitude, it scan still be
hard to identify subtle differences between ESR
spectra, and small imperfections -- such as free MTSL
-- can play an outsized role.

Therefore, here, we use correlation to align the
spectra and then use "dot-product scaling" to normalize
them.  By "dot-product scaling" we mean scaling the
amplitude of one vector (here a spectrum,
:math:`\mathbf{b}`) relative to a reference (here
:math:`\mathbf{c}` to minimize the residual between the
two -- *i.e.* minimizing the following expression by varying the
scaling constant :math:`c`:
:math:`|\mathbf{a}-c\mathbf{b}|^2`
The solution to this is
:math:`c = \frac{\Re[\mathbf{a}\cdot \mathbf{b}]}{||\mathbf{b}||^2}`

In order to do all this, we need a common *x*-axis that
we can use for correlation, etc.
Here, we use the first spectrum in order to determine
the *x*-axis.
Also, for the purposes of dot-product scaling, 
it is better to scale the less noisy spectrum
(:math:`\mathbf{b}` above) 
relative to the noisier spectrum
(:math:`\mathbf{a}` above),
so we assume that the first spectrum is the least
noisy.
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import QESR_scalefactor
from collections import OrderedDict
import matplotlib as mpl
import pickle
init_logging(level='debug')
filenames_w_labels =  [
        ('220307_S175_KCl.DSC','220307_S175_KCl'),
        ('220729_prS175.DSC', '220729 prS175'),
        ('220307_S175_KI.DSC','220307_S175_KI'),                
        ('220307_prS175_KH2PO4.DSC','220307_S175_KH2PO4'), 
        ]
Bname = "$B_0$"
all_files = OrderedDict()
# {{{ load all files first, so that we can determine
# the axis we will want to interpolate onto
for j, (filename, label_str) in enumerate(filenames_w_labels):
    # {{{ load, rescale
    d = find_file(filename, exp_type="francklab_esr/Farhana")
    print(d.name()); quit()
    d /= QESR_scalefactor(d)
    if "harmonic" in d.dimlabels:
        d = d["harmonic", 0]
    d -= d[Bname, :50].C.mean(Bname).data
    if j == 0:
        minB,maxB = d.getaxis(Bname)[r_[0,-1]]
        dB = diff(d.getaxis(Bname)[r_[0,1]]).item()
    else:
        temp = d.getaxis(Bname)[0]
        minB = temp if temp<minB else minB
        temp = d.getaxis(Bname)[-1]
        maxB = temp if temp>maxB else maxB
        temp = diff(d.getaxis(Bname)[r_[0,1]]).item()
        dB = temp if temp<dB else dB
    # }}}
ref_axis = r_[minB:maxB+dB:dB]
# }}}

with figlist_var(width=0.7, filename="ESR_align_example.pdf") as fl:
    # {{{ arrange the figures in the PDF
    fl.par_break()  # each fig on new line
    fl.next("Raw")
    fl.par_break()  # each fig on new line
    fl.next("correlation", legend=True)
    fl.par_break()
    fl.next("aligned, autoscaled", legend=True)
    # }}}
    for label_str, d in all_files:
        # {{{ just show the raw data
        fl.next("Raw")
        fl.plot(d, label=label_str, alpha=0.5)
        # }}}
        if j == 0:
            # if the reference spectrum, store its range
            B_slice = d.getaxis(Bname)[r_[0, -1]]
        d = d.interp(Bname, ref_axis.copy(), kind="linear")
        d[Bname] -= B_start
        if j == 0:
            # u_slice = (-3.41333333,3.41) # hard-coded based on smallest res
            ref_spec_Bdom = d.C
            ref_spec = d.C
            ref_spec.ft(Bname, shift=True)
            ref_spec.run(conj)
            scaling = 1
        else:
            normfactor = d.data.max()
            d.ft(Bname, shift=True)
            correlation = d * ref_spec
            correlation /= normfactor
            correlation.ft_clear_startpoints(
                Bname, t="reset"
            )  # because I want to calculate things in terms of an offset, which can be positive or negative, and need to shift in the next step
            correlation.ift(Bname, shift=True)
            fl.next("correlation")
            fl.plot(correlation, label=label_str)
            thisshift = correlation.real.argmax(Bname).item()
            d *= exp(1j * 2 * pi * d.fromaxis(Bname) * thisshift)
            d.ift(Bname)
            d.set_units(Bname, "G")
            scaling = (d * ref_spec_Bdom).sum(Bname) / (
                ref_spec_Bdom * ref_spec_Bdom
            ).sum(Bname)
            scaling = scaling.real.item()
        d[Bname] += B_start
        d = d[Bname:B_slice]
        fl.next("aligned, autoscaled")
        fl.plot(d / scaling, label=f"{label_str}\nscaling {scaling}", alpha=0.5)
