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
            "9-24-24 13.5 mM TEMPOL sweep",
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
        label = label_str
        fl.basename = "(%s)" % label
        fl.next("raw data")
        fl.image(s, interpolation="auto")
        ###}}}
        ###{{{ Time Domain
        field_axis = []
        for indir_val in range(np.shape(s.getaxis("indirect"))[0]):
            field_axis.append(s.getaxis("indirect")[indir_val][0])
        s.setaxis("indirect", field_axis)
        nu_B12 = s.get_prop("acq_params")["uw_dip_center_GHz"]
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
        fl.image(s["t2":(None, 0.05)])
        ###}}}
        ###{{{ Before Phasing
        # print(np.shape(s['indirect','data{0}','carrierFreq']))
        s = s["t2":(-800, 800)]
        fl.next("before phasing")
        print(psd.ndshape(s))
        if "nScans" in s.dimlabels:
            fl.plot(prscr.select_pathway(s.C.mean('nScans'),signal_pathway))
        else:
            fl.plot(prscr.select_pathway(s, signal_pathway))
        ###}}}
        # fl.show();quit()
        # s.ift('t2')
        # print(s.getaxis('t2')[0])
        # fl.show();quit()
        # s['t2'] -= s.getaxis('t2')[0]
        s.set_units("t2","s")
        s.ft("t2")
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

