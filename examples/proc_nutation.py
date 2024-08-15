"""
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
fig, ax_list = plt.subplots(2, 2)
with psd.figlist_var() as fl:
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle("Nutation %s" % sys.argv[2])
    # }}}
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    # {{{ Apply overall zeroth order correction
    s.ift("t2")
    s /= prscr.zeroth_order_ph(
        prscr.select_pathway(s["t2":0], s.get_prop("coherence_pathway"))
    )
    s.ft("t2")
    fl.image(
        prscr.select_pathway(
            s["t2":signal_range], s.get_prop("coherence_pathway")
        ),
        ax=ax_list[0, 0],
    )
    d_raw = prscr.select_pathway(
        s["t2":signal_range].C, s.get_prop("coherence_pathway")
    )
    ax_list[0, 0].set_title("Signal pathway / ph0")
    # }}}
    # {{{ look at phase variation
    phase_var_s = prscr.select_pathway(
        s["t2":signal_range].C, s.get_prop("coherence_pathway")
    )
    d_uncorrected = phase_var_s.C.real.integrate("t2")
    for j in range(len(s.getaxis("beta"))):
        ph0 = prscr.zeroth_order_ph(phase_var_s["beta", j])
        phase_var_s["beta", j] /= ph0
    d_ind_ph0 = phase_var_s.real.integrate("t2")
    mysign = (d_ind_ph0 / d_uncorrected).angle / np.pi
    mysign = np.exp(1j * np.pi * mysign.run(np.round))
    d_raw *= mysign
    fl.image(d_raw, ax=ax_list[0, 1])
    ax_list[0, 1].set_title("check phase variation along indirect")
    # }}}
    # {{{ apply phasing and FID slice
    s *= mysign
    s.set_error(None)
    s = prscr.fid_from_echo(s, s.get_prop("coherence_pathway"))
    s *= mysign
    s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
    fl.image(s["t2":signal_range], ax=ax_list[1, 0])
    ax_list[1, 0].set_title("Phased and FID sliced")
    # }}}
    s = s["t2":signal_range].real.integrate("t2")
    s.set_error(None)
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real=True)
    fl.plot(s, "o", ax=ax_list[1, 1], human_units=False)
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
    fl.plot(fit, ax=ax_list[1, 1], human_units=False)
    ax_list[1, 1].set_xlabel(r"$\beta$ / $\mathrm{s \sqrt{W}}$")
    ax_list[1, 1].set_ylabel(None)
    beta_90 = f.output("beta_ninety")
    ax_list[1, 1].axvline(beta_90)
    ax_list[1, 1].text(
        beta_90 + 5e-6,
        5e4,
        r"$\beta_{90} = %f \mathrm{\mu s \sqrt{W}}$" % (beta_90 * 1e6),
    )
    ax_list[1, 1].grid()
