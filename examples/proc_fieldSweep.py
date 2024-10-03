"""
Check NMR/ESR resonance ratio using a field sweep
====================================================
Analyzes field sweep data. Determines the optimal field across a gradient that
is on-resonance with the Bridge 12 μw frequency stored in the file to determine
the resonance ratio of MHz/GHz.
"""

import pyspecdata as psd
import pyspecProcScripts as prscr
import os, h5py
import numpy as np
from Instruments.logobj import logobj
import logging
import matplotlib.pyplot as plt

data_target = os.path.normpath(psd.getDATADIR('WK_processed_data'))
signal_pathway = {"ph1": 1}
# {{{ file identifiers 
my_filename = "240924_13p5mM_TEMPOL_field.h5"
my_exptype = "ODNP_NMR_comp/field_dependent"
my_expnode = "field_1"
my_postproctype = "field_sweep_v5"
find_my_file = psd.search_filename(my_filename, exp_type=my_exptype, unique=True)
# }}}
with psd.figlist_var() as fl:
    thisfile, exptype, nodename, post_proc, label_str = (
        my_filename,
        my_exptype,
        my_expnode,
        my_postproctype,
        "240924 13.5 mM TEMPOL field sweep",
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        postproc=post_proc,
        lookup=prscr.lookup_table,
    )
    nu_B12 = s.get_prop("acq_params")["uw_dip_center_GHz"]
    use_freq = True
    if use_freq:
        s["indirect"] = s["indirect"]["carrierFreq"]
        s.set_units("indirect", "MHz")
        s["indirect"] = s["indirect"] / nu_B12
    else:
        # I wanted to use the carrier (use_frep=True), but it seems like the
        # first point isn't stored properly
        s["indirect"] = s["indirect"]["Field"]
        s.set_units("indirect", "G")
    prscr.rough_table_of_integrals(s, fl=fl)
    print(s["indirect"])
    ppt_axis = s["indirect"]
    fl.next("integrated - ppt")
    fl.plot(s.setaxis("indirect", ppt_axis), "o-")
    s.setaxis("indirect", ppt_axis)
    print("s is ", s)
    fitting = s.polyfit("indirect", 4)
    x_min = s["indirect"][0]
    x_max = s["indirect"][-1]
    Field = psd.nddata(np.r_[x_min:x_max:100j], "field")
    fl.plot(Field.eval_poly(fitting, "field"), label="fit")
    print("ESR frequency is %f" % (nu_B12))
    print(
            "The fit finds a max with ppt value:",
            Field.eval_poly(fitting, "field").argmax().item(),
    )    
    print("The data finds a ppt value", abs(s).argmax().item())

