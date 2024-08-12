import pyspecdata as psd
import pyspecProcScripts as prscr
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

signal_range = (-25, 250)
signal_range = (-250, -50)

with psd.figlist_var() as fl:
    s = psd.find_file(
        "240805_amp0p1_27mM_TEMPOL_FID_nutation.h5",
        exp_type="ODNP_NMR_comp/nutation",
        expno="FID_nutation_1",
        lookup=prscr.lookup_table,
    )
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fl.next("Raw data", fig=fig)
    fig.suptitle("Single Pulse nutation")
    signal_pathway = s.get_prop("coherence_pathway")
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    mysgn = prscr.determine_sign(
        prscr.select_pathway(s["t2":signal_range], signal_pathway)
    )
    # {{{ apply overall zeroth order correction
    s /= prscr.zeroth_order_ph(prscr.select_pathway(s["t2":0], signal_pathway))
    d_raw = s.C
    fl.image(prscr.select_pathway(s["t2":signal_range].C, signal_pathway), ax=ax1)
    ax1.set_title("Raw Data")
    d_unc = prscr.select_pathway(s["t2":signal_range].C.real, signal_pathway).integrate(
        "t2"
    )
    # }}}
    # {{{ zeroth_order on individual betas
    for j in range(len(s.getaxis("beta"))):
        ph0 = prscr.zeroth_order_ph(
            prscr.select_pathway(s["beta", j]["t2":signal_range], signal_pathway)
        )
        s["beta", j] /= ph0
    d0 = prscr.select_pathway(s["t2":signal_range].C.real, signal_pathway).integrate(
        "t2"
    )
    # }}}
    # {{{ define mysgn based on phase difference between d_unc and d0
    mysign = (d0 / d_unc).angle / np.pi
    mysign = np.exp(1j * np.pi * mysign.run(np.round))
    d_raw *= mysign
    fl.image(prscr.select_pathway(d_raw, signal_pathway), ax=ax2)
    ax2.set_title("Same Sign")
    s = d_raw
    s.set_units("t2", "s")
    s.set_error(None)
    s = prscr.fid_from_echo(s, signal_pathway)
    s *= mysign
    fl.image(s, human_units=False, ax=ax3)
    ax3.set_title("Phased and FID sliced")
    s = prscr.select_pathway(s["t2":signal_range].real, signal_pathway)
    s.integrate("t2")
    s.set_error(None)
    A, beta_ninety, beta = sp.symbols("A beta_ninety beta", real=True)
    fl.next("Integrated")
    fl.plot(s, "o")
    f = psd.lmfitdata(s)
    f.functional_form = A * sp.sin(beta / beta_ninety * sp.pi / 2)
    f.set_guess(
        A=dict(
            value=s.data.max() * 1.2, min=s.data.max() * 0.8, max=s.data.max() * 1.5
        ),
        beta_ninety=dict(value=2e-5, min=0, max=1),
    )
    f.fit()
    fit = f.eval(100)
    fl.plot(fit)
    psd.gridandtick(plt.gca())
