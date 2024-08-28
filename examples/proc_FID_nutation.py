"""
Process FID nutation data
====================
`py proc_FID_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result as well as examines the phase variation along the indirect dimension.
Finally the data is integrated and fit to a sin function to find the optimal beta_ninety.

Tested with:

``py proc_FID_nutation.py FID_nutation_1 240805_amp0p1_27mM_TEMPOL_FID_nutation.h5 ODNP_NMR_comp/nutation``
"""
import pyspecdata as psd
import pyspecProcScripts as prscr
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
with psd.figlist_var() as fl:
    if "nScans" in s.dimlabels:
        s.mean("nScans")
    s.rename("beta", r"$\beta$")  # label beta with latex
    s, ax_last = prscr.rough_table_of_integrals(
        s, signal_range, fl=fl, title=sys.argv[2], echo_like=False
    )
    # {{{ generate the table of integrals and fit
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety $\\beta$", real=True)
    s = psd.lmfitdata(s)
    s.functional_form = (
        A * sp.exp(-R * beta) * sp.sin(beta / beta_ninety * sp.pi / 2)
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
    fit = s.eval(500)
    fl.plot(fit, ax=ax_last)
    ax_last.set_title("Integrated and fit")
    beta_90 = s.output("beta_ninety")
    ax_last.axvline(beta_90 / 1e-6)
    ax_last.text(
        beta_90 / 1e-6 + 5,
        5e4,
        r"$\beta_{90} = %f\ \mathrm{μs \sqrt{W}}$" % (beta_90 / 1e-6),
    )
    ax_last.grid()
    # }}}