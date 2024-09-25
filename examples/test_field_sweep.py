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
        peakrange = s.get_prop("peakrange")
        print("peakrange is ", peakrange)
        s = s["t2":(-800, 800)]
        fl.next("before phasing")
        print(psd.ndshape(s))
        if "nScans" in s.dimlabels:
            fl.plot(prscr.select_pathway(s.C.mean('nScans'),signal_pathway))
        ###}}}
        # fl.show();quit()
        # s.ift('t2')
        # print(s.getaxis('t2')[0])
        # fl.show();quit()
        # s['t2'] -= s.getaxis('t2')[0]
        s.set_units("t2","s")
        s.ft("t2")
        s = prscr.rough_table_of_integrals(s, fl=fl, title="13.5 mM TEMPOL")
        # fl.show();quit()
        # s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
        fl.next("phased")
        fl.image(s, human_units=False)
        fl.show();quit()
        s = prscr.select_pathway(s, signal_pathway)
        all_offsets = np.zeros(s.shape["indirect"])
        carrier_freq_MHz = []
        for z in range(len(s.getaxis("indirect"))):
            fl.next("Field slicing")
            if "nScans" in s.dimlabels:
                fl.plot(s["indirect", z].C.mean('nScans'), label="scan %d" % z)
                offset = s["indirect", z].C.mean("nScans").argmax("t2").item()
            else:
                offset = s["indirect", z].C.argmax("t2").item()
            all_offsets[z] = offset
            carrier_freq_MHz = s_indirect[:][z]
        s.ift("t2")
        s *= np.exp(
            -1j
            * 2
            * np.pi
            * psd.nddata(all_offsets, [-1], ["indirect"])
            * s.fromaxis("t2")
        )
        s.ft("t2")
        # d = prscr.fid_from_echo(s, signal_pathway, fl=fl)
        peak_range = nu_B12 #s.get_prop("mw_freqs")
        peakrange_prop = s.get_prop("peakrange")
        print("peak_range is", peak_range)
        print("peakrange_prop is", peakrange_prop)
        peak_range = 1.2 * np.r_[-0.5, 0.5] * peak_range
        + np.mean(
            peak_range
        )  # expand by 20%
        frq_slice = peak_range  # s.C.mean('indirect').contiguous(lambda x: x.real > 0.025*s.real.data.max())[0] #dont do, pull prop from FID_from_echo
        # fl.show();quit()
        #frq_slice = peakrange_prop
        fl.next("phased data")
        fl.plot(s)
        plt.axvline(x=frq_slice[0])
        plt.axvline(x=frq_slice[-1])
        if "nScans" in s.dimlabels:
            s = s["t2":frq_slice].C.mean("nScans").integrate("t2")
        else:
            s = s["t2":frq_slice].integrate("t2")
        field_axis_array = np.asarray(field_axis)
        ppt = field_axis_array * s.get_prop("acq_params")["gamma_eff_MHz_G"]
        ppt_axis = np.asarray(ppt)
        ppt_axis /= s.get_prop("acq_params")["uw_dip_center_GHz"]
        # s *= -1
        fl.next("integrated - ppt")
        fl.plot(s.setaxis("indirect", ppt_axis), "o-")
        s.setaxis("indirect", ppt_axis)
        print(s)
        fitting = s.polyfit("indirect", order=4)
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
fl.show()
