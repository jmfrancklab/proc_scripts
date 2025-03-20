"""
Process nutation data
=====================

`py proc_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are
performed according to the `postproc_type` attribute of the data node.
This script plots the result as well as examines the phase variation along the
indirect dimension.
Finally, the data is integrated and fit to a sin**3 function to find the
optimal beta_ninety.

Tested with:

``py proc_nutation.py nutation_1 240805_amp0p1_27mM_TEMPOL_nutation.h5\
        ODNP_NMR_comp/nutation``
"""

import pyspecdata as psd
import pyspecProcScripts as prscr
import sympy as sp
import sys, os
from numpy import r_

# even though the following comes from pint, import the instance from
# pyspecdata, because we might want to mess with it.
from pyspecdata import Q_

if (
    "SPHINX_GALLERY_RUNNING" in os.environ
    and os.environ["SPHINX_GALLERY_RUNNING"] == "True"
):
    sys.argv = [
        sys.argv[0],
        "nutation_1",
        "240805_amp0p1_27mM_TEMPOL_nutation.h5",
        "ODNP_NMR_comp/nutation",
    ]

slice_expansion = 5
assert len(sys.argv) == 4
s = psd.find_file(
    sys.argv[2],
    exp_type=sys.argv[3],
    expno=sys.argv[1],
    lookup=prscr.lookup_table,
)
with psd.figlist_var() as fl:
    frq_center, frq_half = prscr.find_peakrange(prscr.select_pathway(s, s.get_prop("coherence_pathway")), fl=fl)
    signal_range = tuple(slice_expansion * r_[-1, 1] * frq_half + frq_center)
