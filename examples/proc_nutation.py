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
fig, (ax1, ax2) = plt.subplots(1, 2)
with psd.figlist_var() as fl:
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle("Nutation %s" % sys.argv[2])
    # }}}
    signal_pathway = s.get_prop("coherence_pathway")
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    fl.image(
        prscr.select_pathway(s["t2":signal_range], signal_pathway), ax=ax1
    )
    ax1.set_title("Raw Data")
    # {{{ apply overall zeroth order correction
    s /= prscr.zeroth_order_ph(prscr.select_pathway(s["t2":0], signal_pathway))
    s = prscr.fid_from_echo(s, signal_pathway)
    s = prscr.select_pathway(s, signal_pathway)
    fl.image(s["t2":signal_range], human_units=False, ax=ax2)
    ax2.set_title("Phased and FID sliced")
    s = s["t2":signal_range].real.integrate("t2")
    s.set_error(None)
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real=True)
    fl.next("Integrated")
    fl.plot(s, "o")
    f = psd.lmfitdata(s)
    f.functional_form = (
        A * sp.exp(-R * beta) * sp.sin(beta / beta_ninety * sp.pi / 2) ** 3
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
    plt.xlabel(r"$\beta$ / $\mathrm{\mu s \sqrt{W}}$")
    beta_90 = fit.argmax("beta").item() * 1e6
    plt.axvline(beta_90)
    plt.text(
        beta_90 + 5,
        5e4,
        r"$\beta_{90} = %f \mathrm{\mu s \sqrt{W}}$" % beta_90,
    )
    psd.gridandtick(plt.gca())
