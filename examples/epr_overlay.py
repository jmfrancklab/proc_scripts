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
import matplotlib.pylab as plt

mpl.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})

# sphinx_gallery_thumbnail_number = 1

# %% so we can control directories, etc, load the data, but don't mess with it
#    at all (that's handled by align_esr)
filename_dict = {
    "220307_S175_KCl": "220307_S175_KCl.DSC",
    "220729 prS175": "220729_prS175.DSC",
    "220307_S175_KI": "220307_S175_KI.DSC",
    "220307_S175_KH2PO4": "220307_prS175_KH2PO4.DSC",
}
data_dict_multiexpr = {}
for k, v in filename_dict.items():
    data_dict_multiexpr[k] = psd.find_file(v, exp_type="francklab_esr/Farhana")

# %% This one is straightforward, so I just save the output, which I plot
#    below.  This compares several samples with slightly different
#    conditions to look for differences.  Not, in particular, how things
#    get lined up nicely despite the presence of MTSL contaminant in
#    some.
data_dict_multiexpr = align_esr(data_dict_multiexpr)

# %% Next, I load the data from a desalting run (separating protein from
#    MTSL).  I run this at the end.

filename_dict = {}
for j in range(3, 6):
    filename_dict[f"fraction {j}"] = (
        f"240404_L56_MTSL_Rasbatch240320_fraction{j}.DSC"
    )
data_dict_desalt = {}
gen_pdf = False
for k, v in filename_dict.items():
    data_dict_desalt[k] = psd.find_file(v, exp_type="francklab_esr/warren")
# }}}
with psd.figlist_var(width=0.7) as fl:
    fl.next("PR comparison", legend=True, figsize=(3 * 1.05 * 1.618, 3))
    for k, v in data_dict_multiexpr.items():
        fl.plot(v, label=f"{k}\n×{v.get_prop('scaling'):#0.2g}")
    fl.adjust_spines("bottom")
    plt.title("")
    plt.ylabel("")
    plt.gca().set_yticks([])
    align_esr(
        data_dict_desalt,
        fl=fl,
        on_crossing=True,
        correlation_slice=(-0.5e-3, 0.5e-3),
    )
    if gen_pdf:
        fl.show_prep()
        fl.next("centered spectra")
        plt.savefig("overlay.pdf")
