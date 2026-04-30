"""Compare Hall probe field drift with NMR peak shift over time
of a combined_ODNP experiment from log."""

import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt
import re


with psd.figlist_var() as fl:
    for filename, exp_type in [
        ("260406_hydroxytempo_ODNP_1.h5", "B27/ODNP"),  # Broken hdf example
        ("260429_hydroxytempo_ODNP_2.h5", "B27/ODNP"),
    ]:
        s = psd.find_file(
            filename,
            exp_type=exp_type,
            expno="ODNP",
            lookup=prscr.lookup_table,
        )
        frq_center, frq_half = prscr.find_peakrange(
            s, direct="t2", peak_lower_thresh=0.05
        )
        frq_range = (
            frq_center - frq_half,
            frq_center + frq_half,
        )
        if s.get_prop("log") is None:
            s = prscr.attach_log_data_from_file(s, filename, exp_type)
        # {{{ gen coords if old data
        m = re.search(".*ODNP.*v([0-9]+)$", s.get_prop("postproc_type"))
        if m:
            vernum = int(m.groups()[0])
        else:
            raise IOError(
                "What the heck type of postproc type is that!!"
                f" ({s.get_prop('postproc_type')})"
            )
        if vernum < 6:
            # the next line is done already if we are 6 or higher
            s = prscr.generate_coordinates_from_log.generate_coordinates_from_log(
                s
            )
        # }}}
        log_array = s.get_prop("log").total_log
        log_start_time = log_array["time"][0].item()
        log_array["time"] -= log_start_time
        Hall_drift_Hz = (
            psd.nddata(
                1e6
                * (
                    log_array["field"]
                    * s.get_prop("acq_params")["gamma_eff_MHz_G"]
                    - s.get_prop("acq_params")["carrierFreq_MHz"]
                ),
                [-1],
                ["t"],
            )
            .setaxis("t", log_array["time"])
            .set_units("t", "s")
        )
        s["indirect"] = s["indirect"]["time"]
        s.set_error("indirect", 0)
        s.rename("indirect", "t").set_units("t", "s")
        # the time axis for the signal is originally a structured array that gives
        # both experiment start and stop times.  Convert to a normal array
        # {{{ rather than averaging over nScans, I create an appropriate time
        #     axis for it, and then smoosh to create a dimension with all my scans
        if "nScans" in s.dimlabels:
            s["nScans"] = s["nScans"] * (
                s.get_prop("acq_params")["acq_time_ms"] * 1e-3
                + s.get_prop("acq_params")["repetition_us"] * 1e-6
            )
            s.smoosh(["t", "nScans"], dimname="t")
            # smoosh makes a structured array, which I convert to a normal array:
            s["t"] = s["t"]["nScans"] + s["t"]["t"]
        s.set_units("t", "s")
        # }}}
        s.reorder("t2", first=False)
        fl.basename = filename
        fig = fl.next("NMR signal - $\\varphi_0$ only")
        fig.clf()
        s /= prscr.zeroth_order_ph(s)
        fl.DCCT(s)
        fl.next("slice FID")
        s = prscr.fid_from_echo(
            s, signal_pathway=s.get_prop("coherence_pathway")
        )
        s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
        fl.plot(s.real)
        fl.next("NMR signal - with zf and conv (tdom)")
        # the linewidth is about 1 kHz, so in the following, we're shooting
        # for a matched filter
        s.ift("t2").ft("t2", pad=s.shape["t2"] * 20).convolve("t2", 1e3)
        fl.plot(s)
        fl.next("normalized (zoomed)")
        s = s["t2":frq_range].real
        s /= s.C.integrate("t2")
        fl.plot(s)
        # {{{ now that we have s(ν)/∫s(ν')dν', calculate
        #     ∫ ν (s(ν)/∫s(ν')dν') dν
        s *= s.fromaxis("t2")
        s.integrate("t2")
        # }}}
        fl.next("B field and peak shift vs. time")
        fl.plot(Hall_drift_Hz, ".", label="Hall probe")
        fl.plot(s, "o", label="NMR")
        plt.ylim(frq_range)
