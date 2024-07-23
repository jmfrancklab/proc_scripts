"""
Process nutation data
====================
`py proc_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.

Tested with:

``py proc_nutation.py nutation_2 240710_27mM_TEMPOL_chokes_T_atprobe_nutation.h5 ODNP_NMR_comp/nutation``

"""
import pyspecdata as psd
import pyspecProcScripts as prscr
import numpy as np
import matplotlib.pyplot as plt
import sys

assert len(sys.argv) == 4
d = psd.find_file(
    sys.argv[2],
    exp_type=sys.argv[3],
    expno=sys.argv[1],
    lookup=prscr.load_data.lookup_table,
)
if d.get_prop("postproc_type") == "spincore_SE_v1":
    d.rename("indirect", "p_90")
    d["p_90"] *= 1e-9
elif "indirect" in d.dimlabels:
    d.rename("indirect", "p_90")
signal_range = (-100, 100)
disprange = (-500, 500)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
with psd.figlist_var() as fl:
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle("Nutation %s" % sys.argv[2])
    # }}}
    # {{{ average scans and show raw
    if "nScans" in d.dimlabels:
        d.mean("nScans")
    d = d["t2":disprange]
    og_signs = prscr.determine_sign(prscr.select_pathway(d["t2":signal_range]))
    d.reorder("t2", first=False)
    fl.image(d.C["t2":disprange], interpolation="auto", ax=ax1)
    ax1.set_title("Raw Data")
    # }}}
    # {{{ apply clock correction
    total_corr = 0
    for j in range(5):
        corr = prscr.clock_correction(
            prscr.select_pathway(d, d.get_prop("coherence_pathway"))
            * og_signs
            * np.exp(-1j * 2 * np.pi * total_corr * d.fromaxis("p_90")),
            "p_90",
        )
        total_corr += corr
    d *= np.exp(-1j * 2 * np.pi * total_corr * d.fromaxis("p_90"))
    # }}}
    # {{{ apply zeroth order correction for final sign assignment
    for j in range(len(d.getaxis("p_90"))):
        ph0 = prscr.zeroth_order_ph(
            prscr.select_pathway(d["p_90", j]["t2":signal_range])
        )
        d["p_90", j] /= ph0
    my_signs = prscr.determine_sign(prscr.select_pathway(d["t2":signal_range]))
    # }}}
    d = prscr.fid_from_echo(d, d.get_prop("coherence_pathway"))
    align_option = input("Do you want to apply the correlation alignment?")
    if align_option.lower().startswith("n"):
        pass
    elif align_option.lower().startswith("y"):
        d.ift("t2")
        lambda_L = prscr.fit_envelope(d)
        d.ft("t2")
        align_sign = prscr.determine_sign(
            prscr.select_pathway(d, d.get_prop("coherence_pathway"))
        )
        matched = (prscr.select_pathway(d) * align_sign).ift("t2")
        matched *= np.exp(-np.pi * lambda_L * matched.fromaxis("t2"))
        matched.ft("t2")
        frq_atmax = matched.real.argmax("t2")
        d.ift("t2")
        d *= np.exp(-1j * 2 * np.pi * frq_atmax * d.fromaxis("t2"))
        d.ft("t2")
        d.ift(list(d.get_prop("coherence_pathway")))
        opt_shift, sigm, mask = prscr.correl_align(
            d,
            indirect_dim="p_90",
            sigma=lambda_L,
            signal_pathway=d.get_prop("coherence_pathway"),
        )
        d.ift("t2")
        d *= np.exp(-1j * 2 * np.pi * opt_shift * d.fromaxis("t2"))
        d.ft(list(d.get_prop("coherence_pathway")))
        d = d["t2":(0, None)]
        d *= 2
        d["t2":0] *= 0.5
        d.ft("t2")
    fl.image(d, interpolation="auto", ax=ax2)
    ax2.set_title("Phased and FID sliced")
    d *= og_signs
    d.set_units("p_90", "s")
    fl.image(d, interpolation="auto", ax=ax3)
    ax3.set_title("Flip sign back")
