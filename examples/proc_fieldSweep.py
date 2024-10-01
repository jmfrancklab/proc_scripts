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
from Instruments.logobj import logobj
import logging
data_target = os.path.normpath(psd.getDATADIR('WK_processed_data'))

signal_pathway = {"ph1": 1}
# {{{ file identifiers 
my_filename = "240924_13p5mM_TEMPOL_field.h5"
my_exptype = "ODNP_NMR_comp/field_dependent"
my_expnode = "field_1"
my_postproctype = "field_sweep_v4"
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
    use_freq = False
    if use_freq:
        s["indirect"] = s["indirect"]["carrierFreq"] * 1e6
        s.set_units("indirect", "Hz")
    else:
        # I wanted to use the carrier (use_frep=True), but it seems like the
        # first point isn't stored properly
        s["indirect"] = s["indirect"]["Field"]
        s.set_units("indirect", "G")
    prscr.rough_table_of_integrals(s, fl=fl)
    # {{{Obtain ESR frequency and chunk/reorder dimensions	
    nu_B12_GHz = (s.get_prop("acq_params")["mw_freqs"][0]) / 1e9	
    print(nu_B12_GHz)
    # }}}	
    # {{{DC offset correction	
    s.ift("t2")	
    s.ift("ph1")	
    t2_max = s.getaxis("t2")[-1]	
    rx_offset_corr = s["t2" : (t2_max * 0.75, None)]	
    rx_offset_corr = rx_offset_corr.mean(["t2"])	
    s -= rx_offset_corr	
    s.ft(["ph1"])	
    s.ft("t2")	
    # }}}	
    # {{{ set up figure and plot raw data	
    fig, ax_list = subplots(1, 4, figsize=(10, 3))	
    fl.next("Field Sweep Processing", fig=fig)	
    fl.image(s, ax=ax_list[0])	
    ax_list[0].set_title("Raw data\nFrequency Domain")	
    # }}}	
    # {{{frequency filtering and phase correct	
    s = s["t2":(-1e3, 1e3)]	
    s.ift("t2")	
    best_shift = hermitian_function_test(select_pathway(s, signal_pathway))	
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})	
    s = s["t2":(0, None)]	
    s["t2", 0] *= 0.5	
    s.ft("t2")	
    s = select_pathway(s, signal_pathway)	
    s /= zeroth_order_ph(s.C.mean("t2"))	
    nu_NMR = []	
    assert set(s.getaxis("indirect").dtype.names) == {	
        "Field",	
        "carrierFreq",	
    }, "'indirect' axis should be a structured array that stores the carrier frequency and the field"	 

