""" Calculate covariance matrix for NMR signal
==============================================
Here, we show that by applying a correlation alignment we drastically decrease
the covariance of NMR signal. This is due to the fact that the peak is moving
around in frequency space -- the phase of the time domain points is not fixed.
As the phase varies, the real and imaginary components of the signal 
*change together* - covariance.
"""
import pyspecdata as psd
import numpy as np
import pyspecProcScripts as psdpr
import matplotlib.pyplot as plt

signal_pathway = {"ph1": 1}
d = psd.find_file(
    "241003_27mM_TEMPOL_amp0p1_var_tau_pm_echo",
    exp_type="ODNP_NMR_comp/Echoes",
    expno="echo_1",
    postproc="none",
).squeeze()
# {{{ Make copy of data and align to compare with unaligned
#     data
aligned = d.C
aligned.ft("t2", shift=True)
opt_shift, sigma, mask_fn = psdpr.correl_align(
    aligned,
    indirect_dim="nScans",
    signal_pathway=signal_pathway,
)
aligned.ift("t2")
aligned *= np.exp(-1j * 2 * np.pi * opt_shift * aligned.fromaxis("t2"))
# }}}
with psd.figlist_var() as fl:
    for thisdata, label in [(d, "Unaligned"), (aligned, "With Alignment")]:
        thisdata = thisdata["ph1", 0]  # pull only one phcyc step to simplify
        # I explictly reorder so I know which dimension is which when I do raw
        # numpy and matplotlib
        thisdata.reorder(["t2", "nScans"])
        fig, ax_list = plt.subplots(3)
        fig.suptitle(label)
        fl.next("%s" % label, fig=fig)
        fl.image(abs(thisdata.cov_mat("nScans")), ax=ax_list[0])
        ax_list[0].set_title("Full covariance")
        # now I zoom in to show that when I'm in the part of the FID that's
        # just noise (no shifting frequencies, I get the diagonal covariance
        # that I expect)
        fl.image(
            abs(thisdata["$t2_i$":(0.6, 0.7)]["$t2_j$":(0.6, 0.7)]),
            ax=ax_list[1],
        )
        ax_list[1].set_title("covariance, zoomed late times")
        fl.image(
            abs(thisdata["$t2_i$":(0.9, 1.0)]["$t2_j$":(0.9, 1.0)]),
            ax=ax_list[2],
        )
        ax_list[2].set_title("covariance, zoomed even later")
