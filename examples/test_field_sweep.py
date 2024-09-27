import pyspecdata as psd
import pyspecProcScripts as prscr
import numpy as np
import matplotlib.pyplot as plt

signal_pathway = {"ph1": 1}
with psd.figlist_var() as fl:
    for thisfile, exptype, nodename, post_proc, label_str in [
        (
            "240924_13p5mM_TEMPOL_field.h5",
            "ODNP_NMR_comp/field_dependent",
            "field_1",
            "field_sweep_v4",
            "240924 13.5 mM TEMPOL field sweep",
        )
    ]:
        s = psd.find_file(
            thisfile,
            exp_type=exptype,
            expno=nodename,
            postproc=post_proc,
            lookup=prscr.lookup_table,
        )
        s_indirect = s["indirect"]
        ###{{{ Raw data
        fl.basename = "(%s)" % label_str
        fl.next("raw data")
        fl.image(s, interpolation="auto")
        ###}}}
        ###{{{ Time Domain
        print(s)
        print(psd.ndshape(s))
        #fl.show();quit()
        field_axis = s["indirect"]        
        s.setaxis("indirect", field_axis)
        print("field_axis is", field_axis)
        nu_B12, half_window = prscr.find_peakrange(s, fl=fl)
        if "nScans" in s.dimlabels:
            s.reorder(["ph1", "indirect", "nScans", "t2"])
        else:
            s.reorder(["ph1", "indirect", "t2"])
        s.ift("t2")
        s.ift(["ph1"])
        t2_max = s.getaxis("t2")[-1]
        rx_offset_corr = s["t2" : (t2_max * 0.75, None)]
        rx_offset_corr = rx_offset_corr.mean(["t2"])
        s -= rx_offset_corr
        fl.next("time domain")
        s.ft(["ph1"])
        fl.image(s["t2":(None, t2_max)])
        print("t2 units", s.get_units("t2"))
        #fl.show();quit()
        ###}}}
        ###{{{ Before Phasing
        #s.ft("t2")
        s = s["t2":s.get_prop("peakrange")]
        fl.next("before phasing")
        if "nScans" in s.dimlabels:
            fl.plot(prscr.select_pathway(s.C.mean('nScans'),signal_pathway))
        print("peakrange before rough func", s.get_prop("peakrange"))
        s.ft("t2")
        s = prscr.rough_table_of_integrals(s, fl=fl, title=label_str)
        print("peakrange after rough func", s.get_prop("peakrange"))
        fl.next("before phasing")
        print(psd.ndshape(s))
        fl.show();quit()
        ###}}}
        s.set_units("t2","s")
        s = prscr.fid_from_echo(s, signal_pathway, fl=fl)
        print(s.get_prop("peakrange"))
        peakrange = s.get_prop("peakrange")
        print("peakrange is ", peakrange)
        print("s shape after fid_from_echo", psd.ndshape(s))
        s, _ = prscr.rough_table_of_integrals(s, fl=fl, title="13.5 mM TEMPOL")
        print(s)
        print("shape of s", psd.ndshape(s))
        #fl.show();quit()
        # s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
        fl.next("phased")
        fl.plot(s, 'gX', human_units=False)
fl.show()

