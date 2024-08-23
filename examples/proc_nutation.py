"""
Process nutation data
====================
`py proc_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result as well as examines the phase variation along the indirect dimension.
Finally, the data is integrated and fit to a sin**3 function to find the optimal beta_ninety.

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
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_figwidth(15)
fig.set_figheight(6)
with psd.figlist_var() as fl:
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    # {{{ set up subplots
    fl.next("Raw Data with averaged scans", fig=fig)
    fig.suptitle("Nutation %s" % sys.argv[2])
    # }}}
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
        ax=ax1,
        human_units=False,
    )
    ax1.set_title("Signal pathway / ph0")
    # }}}
    # {{{ Check phase variation along indirect
    mysign = prscr.determine_sign(
        s, "beta", signal_range
    )
    s *= mysign
    fl.image(
        prscr.select_pathway(
            s["t2":signal_range], s.get_prop("coherence_pathway")
        ),
        ax=ax2,
        human_units=False
    )
    ax2.set_title("Check phase variation along indirect")
    # }}}
    # {{{ apply phasing and FID slice
    s = prscr.fid_from_echo(s.set_error(None), s.get_prop("coherence_pathway"))
    s *= mysign
    s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
    fl.image(s["t2":signal_range], ax=ax3)
    ax3.set_title("Phased and FID sliced")
    # }}}
    # {{{ generate the table of integrals and fit
    s = s["t2":signal_range].real.integrate("t2").set_error(None)
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real=True)
    fl.next("Integrated and Fit")
    fl.plot(s, "o")
    s = psd.lmfitdata(s)
    s.functional_form = (
        A * sp.exp(-R * beta) * sp.sin(beta / beta_ninety * sp.pi / 2) ** 3
    )
    s.set_guess(
        A=dict(
            value=s.data.max() * 1.2,
            min=s.data.max() * 0.8,
            max=s.data.max() * 1.5,
        ),
        R=dict(value=3e3, min=0, max=3e4),
        beta_ninety=dict(value=2e-5, min=0, max=1),
    )
    s.fit()
    fit = s.eval(100)
    fl.plot(fit)
    plt.xlabel(r"$\beta$ / $\mathrm{s \sqrt{W}}$")
    plt.ylabel(None)
    beta_90 = s.output("beta_ninety")
    plt.axvline(beta_90)
    plt.text(
        beta_90 + 5,
        5e4,
        r"$\beta_{90} = %f \mathrm{s \sqrt{W}}$" % beta_90,
    )
    psd.gridandtick(plt.gca())
    # }}}
