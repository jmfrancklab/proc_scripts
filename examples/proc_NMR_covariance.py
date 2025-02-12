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
    "241010_ssProbe_toroid_pmp90_16p7_echo",
    exp_type="ODNP_NMR_comp/Echoes",
    expno="echo_1",
    lookup=psdpr.lookup_table,
)
# {{{ Make copy of data and align to compare with unaligned data
aligned = d.C
aligned.ift("ph1")  # go into phase cycling dimension for alignment
opt_shift, sigma, mask_fn = psdpr.correl_align(
    aligned,
    indirect_dim="nScans",
    signal_pathway=signal_pathway,
)
aligned.ift("t2")
aligned *= np.exp(-1j * 2 * np.pi * opt_shift * aligned.fromaxis("t2"))
aligned.ft("ph1")
# }}}
d.ift("t2")
with psd.figlist_var() as fl:
    for thisdata, label in [(d, "Unaligned"), (aligned, "With Alignment")]:
        thisdata = thisdata["ph1", 0]  # pull only one phcyc step to simplify
        # {{{ Set up figure
        fig, ax_list = plt.subplots(3)
        fig.suptitle(label)
        fl.next("%s" % label, fig=fig)
        ax_list[0].set_title("Covariance-variance matrix of full echo")
        ax_list[1].set_title(
            "Covariance-variance matrix, zoomed into late times"
        )
        ax_list[2].set_title(
            "Covariance-variance matrix, zoomed into even later times"
        )
        # }}}
        # Show cov-var matrix of full time domain
        fl.image(abs(thisdata.cov_mat("nScans")), ax=ax_list[0])
        # now I zoom in (in two steps) to show that when I'm in the part of the
        # FID that's just noise (no shifting frequencies, I get the diagonal
        # covariance that I expect)
        fl.image(
            abs(thisdata["$t2_i$":(0.6, 0.7)]["$t2_j$":(0.6, 0.7)]),
            ax=ax_list[1],
        )
        fl.image(
            abs(thisdata["$t2_i$":(0.9, 1.0)]["$t2_j$":(0.9, 1.0)]),
            ax=ax_list[2],
        )
