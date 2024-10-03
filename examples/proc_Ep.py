"""
Process Enhancement experiment 
====================================================
Opens .h5 results file, uses rough_table_of_integrals() to roughly process
dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
import os

data_target = os.path.normpath(psd.getDATADIR("WK_processed_data"))

with psd.figlist_var() as fl:
    thisfile, exptype, nodename, lookup = (
        "240924_13p5mM_TEMPOL_ODNP_1.h5",
        "ODNP_NMR_comp/ODNP",
        "ODNP",
        prscr.lookup_table,
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        lookup=prscr.lookup_table,
    )
    s["indirect"] = s["indirect"]["start_times"]
    s.set_units("indirect", "s")
    prscr.rough_table_of_integrals(s, fl=fl)
