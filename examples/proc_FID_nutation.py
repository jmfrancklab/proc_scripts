"""
=====================
Process FID nutation data
====================
`py proc_FID_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.

Tested with:

``py proc_FID_nutation.py FID_nutation_1 240805_amp0p1_27mM_TEMPOL_FID_nutation.h5 ODNP_NMR_comp/nutation``
"""
import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
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
fig.set_figwidth(15)
fig.set_figheight(6)
with psd.figlist_var() as fl:
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle("FID Nutation %s" % sys.argv[2])
    # }}}
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    # {{{ apply overall zeroth order correction
    s = prscr.select_pathway(
        s["t2":signal_range], s.get_prop("coherence_pathway")
    )
    s.ift("t2")
    s /= prscr.zeroth_order_ph(s["t2":0])
    s.ft("t2")
    fl.image(s, ax=ax1)
    d_raw = s.C
    ax1.set_title("Signal pathway / ph0")
    # }}}
    # {{{ Look at phase variation
    phase_var_s = s.C
    d_uncorrected = phase_var_s.C.real.integrate("t2")
    for j in range(len(s.getaxis("beta"))):
        ph0 = prscr.zeroth_order_ph(phase_var_s["beta", j])
        phase_var_s["beta", j] /= ph0
    d_ind_ph0 = phase_var_s.real.integrate("t2")
    mysign = (d_ind_ph0 / d_uncorrected).angle / np.pi
    mysign = np.exp(1j * np.pi * mysign.run(np.round))
    d_raw *= mysign
    fl.image(d_raw, ax=ax2)
    ax2.set_title("Check phase variation along indirect")
    # }}}
    s = s.real.integrate("t2")
    s.set_error(None)
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real=True)
    fl.plot(s, "o", ax=ax3, human_units=False)
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
    fl.plot(fit, ax=ax3, human_units=False)
    ax3.set_title("Integrated and fit")
    ax3.set_xlabel(r"$\beta$ / $\mathrm{\mu s \sqrt{W}}$")
    beta_90 = fit.argmax("beta").item() * 1e6
    ax3.axvline(beta_90)
    ax3.text(
        beta_90 + 5,
        5e4,
        r"$\beta_{90} = %f \mathrm{\mu s \sqrt{W}}$" % beta_90,
    )
    ax3.grid()
