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

from pyspecProcScripts import align_esr
import matplotlib as mpl
import pyspecdata as psd

mpl.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})

# sphinx_gallery_thumbnail_number = 4

# {{{ so we can control directories, etc, load the data, but don't mess with it
#     at all (that's handled by align_esr)
filename_dict = {
    "220307_S175_KCl": "220307_S175_KCl.DSC",
    "220729 prS175": "220729_prS175.DSC",
    "220307_S175_KI": "220307_S175_KI.DSC",
    "220307_S175_KH2PO4": "220307_prS175_KH2PO4.DSC",
}
data_dict = {}
for k, v in filename_dict.items():
    data_dict[k] = psd.find_file(v, exp_type="francklab_esr/Farhana")
# }}}
with psd.figlist_var(width=0.7, filename="ESR_align_example.pdf") as fl:
    align_esr(data_dict, fl=fl)
    # fl.next("centered spectra")
    # plt.savefig("single_mutant_overlay.pdf")
