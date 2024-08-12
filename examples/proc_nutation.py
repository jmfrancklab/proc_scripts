"""
=====================
Process nutation data
====================
`py proc_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.

Tested with:

``py proc_nutation.py nutation_1 240805_amp0p1_27mM_TEMPOL_nutation.h5 ODNP_NMR_comp/nutation``
"""
import pyspecdata as psd
import pyspecProcScripts as prscr
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import sys

signal_range = (-250, 250)
assert len(sys.argv) == 4
s = psd.find_file(
    sys.argv[2],
    exp_type=sys.argv[3],
    expno=sys.argv[1],
    lookup=prscr.lookup_table,
)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
with psd.figlist_var() as fl:
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle("Nutation %s" % sys.argv[2])
    # }}}
    signal_pathway = s.get_prop("coherence_pathway")
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    mysgn = prscr.determine_sign(
        prscr.select_pathway(s["t2":signal_range], signal_pathway)
    )
    # {{{ apply overall zeroth order correction
    s /= prscr.zeroth_order_ph(prscr.select_pathway(s["t2":0], signal_pathway))
    d_raw = s.C
    fl.image(
        prscr.select_pathway(s["t2":signal_range].C, signal_pathway), ax=ax1
    )
    ax1.set_title("Raw Data")
    d_unc = prscr.select_pathway(
        s["t2":signal_range].C.real, signal_pathway
    ).integrate("t2")
    # }}}
    # {{{ zeroth_order on individual betas
    for j in range(len(s.getaxis("beta"))):
        ph0 = prscr.zeroth_order_ph(
            prscr.select_pathway(
                s["beta", j]["t2":signal_range], signal_pathway
            )
        )
        s["beta", j] /= ph0
    d0 = prscr.select_pathway(
        s["t2":signal_range].C.real, signal_pathway
    ).integrate("t2")
    # }}}
    # {{{ define mysgn based on phase difference between d_unc and d0
    mysign = (d0 / d_unc).angle / np.pi
    mysign = np.exp(1j * np.pi * mysign.run(np.round))
    # }}} 
    d_raw *= mysign
    fl.image(prscr.select_pathway(d_raw, signal_pathway), ax=ax2)
    ax2.set_title("Same Sign")
    s = d_raw
    s.set_units("t2", "s")
    s.set_error(None)
    s = prscr.fid_from_echo(s, signal_pathway)
    if mysign["beta", 5] < 0:
        s *= -mysign
    else:
        s *= mysign
    fl.image(s, human_units=False, ax=ax3)
    ax3.set_title("Phased and FID sliced")
    s = prscr.select_pathway(s["t2":signal_range].real, signal_pathway)
    s.integrate("t2")
    s.set_error(None)
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real=True)
    fl.next("Integrated")
    fl.plot(s, "o")
    f = psd.lmfitdata(s)
    f.functional_form = (
        A * sp.exp(-R * beta) * sp.sin(beta / beta_ninety * sp.pi / 2)
    )
    f.set_guess(
        A=dict(
            value=s.data.max() * 1.2,
            min=s.data.max() * 0.8,
            max=s.data.max() * 1.5,
        ),
        R=dict(value=3e3, min=0, max=3e4),
        beta_ninety=dict(value=2e-5, min=0, max=1),
    )
    f.fit()
    fit = f.eval(100)
    fl.plot(fit)
    psd.gridandtick(plt.gca())
