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
    
    field_axis = s["indirect"]
    nu_B12_GHz = (s.get_prop("acq_params")["uw_dip_center_GHz"]) / 1e9	
    print("B12 freq is ", nu_B12_GHz)
    s.ift("t2")	
    s.ift("ph1")	
    t2_max = s.getaxis("t2")[-1]	
    rx_offset_corr = s["t2" : (t2_max * 0.75, None)]	
    rx_offset_corr = rx_offset_corr.mean(["t2"])	
    s -= rx_offset_corr	
    s.ft(["ph1"])	
    s.ft("t2")	
    fl.next("Field Sweep Processing")	
    fl.image(s)	
    s = s["t2":(-1e3, 1e3)]	
    s.ift("t2")	
    s = s["t2":(0,None)]
    best_shift = prscr.hermitian_function_test(prscr.select_pathway(s, signal_pathway))	
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})	
    #s = s["t2":(0, None)]	
    s["t2", 0] *= 0.5	
    s.ft("t2")	
    s = prscr.select_pathway(s, signal_pathway)	
    s /= prscr.zeroth_order_ph(s.C.mean("t2"))	
    print("s's dimlabels are ", s.dimlabels)
    #assert set(s.getaxis("indirect").dtype.names) == {	
    #    "Field",	
    #    "carrierFreq",	
    #}, "'indirect' axis should be a structured array that stores the carrier frequency and the field"	 
    all_offsets = np.zeros(len(s.getaxis("indirect")))
    carrier_freq_MHz = []
    for z in range(len(s.getaxis("indirect"))):
        fl.next("Field slicing")
        if "nScans" in s.dimlabels:
            fl.plot(s["indirect", z].C.mean("nScans"), label="scan %d" % z)
            offset = s["indirect", z].C.mean("nScans").argmax("t2").item()
        else:
            fl.plot(s["indirect", z], label="scan %d" % z)
            offset = s["indirect", z].C.argmax("t2").item()
        all_offsets[z] = offset
        carrier_freq_MHz = s["indirect"][:][z]
    s.ift("t2")
    s *= np.exp(
        -1j
        * 2
        * np.pi
        * psd.nddata(all_offsets, [-1], ["indirect"])
        * s.fromaxis("t2")
    )    
    s.ft("t2")
    peak_range = nu_B12_GHz
    peak_range = 1.2 * np.r_[-0.5, 0.5] * peak_range + np.mean(peak_range)
    frq_slice = peak_range
    fl.next("phased data")
    fl.plot(s)
    plt.axvline(x=frq_slice[0])
    plt.axvline(x=frq_slice[-1])
    print("t2 as a dimlabel is ", s["t2"])
    if "nScans" in s.dimlabels:
        s = s["t2":frq_slice].C.mean("nScans").integrate("t2")
    else:
        s = s["t2":frq_slice].integrate("t2")
    field_axis_array = np.asarray(field_axis)
    ppt = field_axis_array * s.get_prop("acq_params")["gamma_eff_MHz_G"]
    ppt_axis = np.asarray(ppt)
    ppt_axis /= s.get_prop("acq_params")["uw_dip_center_GHz"]
    fl.next("integrated - ppt")
    fl.plot(s.setaxis("indirect", ppt_axis), "o-")
    s.setaxis("indirect", ppt_axis)
    print(s)
    fl.show();quit()
    fitting = s.polyfit("indirect", order=4)
    x_min = s["indirect"][0]
    x_max = s["indirect"][-1]
    Field = psd.nddata(np.r_[x_min:x_max:100j], "field")
    fl.plot(Field.eval_poly(fitting, "field"), label="fit")
    print("ESR frequency is %f" % (nu_B12_GHz))
    print(
        "The fit finds a max with ppt value:",
        Field.eval_poly(fitting, "field").argmax().item(),
    )
    print("The data finds a ppt value", abs(s["indirect"]).argmax().item())
   
