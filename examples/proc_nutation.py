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
    s, ax3 = prscr.rough_table_of_integrals(
        s, signal_range, fl=fl, echo_like=True, title=sys.argv[2]
    )
    # {{{ generate the table of integrals and fit
    A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real=True)
    fl.plot(s, "o", ax=ax3)
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
    fit = s.eval(500)
    fl.plot(fit, ax=ax3)
    ax3.set_xlabel(r"$\beta$ / $\mathrm{\mu s \sqrt{W}}$")
    ax3.set_ylabel(None)
    beta_90 = s.output("beta_ninety") / 1e-6
    ax3.axvline(beta_90)
    ax3.text(
        beta_90 + 5,
        5e4,
        r"$\beta_{90} = %f \mathrm{s \sqrt{W}}$" % beta_90,
    )
    ax3.grid()
    # }}}
