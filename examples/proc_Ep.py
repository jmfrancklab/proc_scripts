"""
Check
====================================================
Analyzes
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
from Instruments.logobj import logobj

signal_pathway = {"ph1": 1}
with psd.figlist_var() as fl:
    thisfile, exptype, nodename, post_proc, lookup = (
        "240924_13p5mM_TEMPOL_ODNP_1.h5",
        "ODNP_NMR_comp/ODNP",
        "ODNP",
        "spincore_ODNP_v5",
        prscr.lookup_table,
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        postproc=post_proc,
        lookup=prscr.lookup_table,
    )
    s["indirect"] = s["indirect"]["start_times"]
    s.set_units("indirect", "s")
    prscr.rough_table_of_integrals(s, fl=fl)
    fl.show()

