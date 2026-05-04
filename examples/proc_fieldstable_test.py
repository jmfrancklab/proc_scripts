# TODO ☐: I don't understand the point of this file.  We already have
#         plot_field.  Why can't we just use that to process these types
#         of files?
"""
Process repeated echo field-stability test data
===============================================

# TODO ☐: so, this name was new, so in FLInst, I just said "well,
#         everything that's called stability_test will be newer data, so
#         it should just be v1, not v4.
#         But, then I see this, and don't understand it.  Do you have
#         files with these versions that you actually want? Or are these
#         just things that you acquired over the past week, in which
#         case, I would say, let's just rename the postproc type of
#         those to "Unsupported."  Overall, I'm just confused about what
#         the point of this script is given that plot_fields already
#         exists.
``stability_test_v1`` files do not store enough time/field information for
this script.  ``stability_test_v2`` files have no log, ``stability_test_v3``
files have a log that is not attached, and ``stability_test_v4`` files
already have both the log and start/stop times.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd
import pyspecProcScripts as prscr
import re
from numpy import pi
from pyspecProcScripts.generate_coordinates_from_log import (
    generate_coordinates_from_log,
)


plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
plt.rcParams.update(
    {
        "errorbar.capsize": 2,
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.figsize": (6, 4),
    }
)

with psd.figlist_var() as fl:
    for filename, exp_type, nodename in [
        ("251110_hydroxytempo_n_scan.h5", "B27/n_scans", "n_scan_37"),  # v1
        ("260409_hydroxytempo_n_scan.h5", "B27/n_scans", "n_scan_1"),  # v2
        (
            "260414_hydroxytempo_n_scan.h5",
            "B27/n_scans",
            "n_scan_1",
        ),  # v3 distrupted log
        ("260414_hydroxytempo_n_scan.h5", "B27/n_scans", "n_scan_2"),  # v3
    ]:
        fl.basename = f"{filename} {nodename}"
        s = psd.find_file(
            filename,
            exp_type=exp_type,
            expno=nodename,
            lookup=prscr.lookup_table,
            fl=None,
        )
        postproc_type = s.get_prop("postproc_type")
        print("using postproc type", postproc_type)
        # {{{ gen coords if old data
        m = re.search(".*stability_test.*v([0-9]+)$", postproc_type)
        if m:
            vernum = int(m.groups()[0])
        else:
            raise IOError(
                "What the heck type of postproc type is that!!"
                f" ({postproc_type})"
            )
        if vernum == 1:
            print(
                f"This script does not work for {postproc_type} files, since"
                f" the field and time is not stored at all, so I won't plot "
                f"{filename}"
            )
            continue
        old_axis = s["indirect"].copy()
        if vernum > 2 and s.get_prop("log") is None:
            # the log has not changed, just its location in the HDF file.
            s = prscr.attach_log_data_from_file(s, filename, exp_type)
        if s.get_prop("log") is not None and not hasattr(
            s.get_prop("log"), "total_log"
        ):
            # Convert a stored log to an logobj, if it hasn't.
            s.set_prop("log", prscr.logobj.from_group(s.get_prop("log")))
        use_fake_log = vernum == 2
        if vernum == 3:
            log_times = s.get_prop("log").total_log["time"]
            # If the log times do not cover the old axis times, then we need to
            # make a fake log. This was due to sava data does not write the
            # experiment number to the file name.
            use_fake_log = (
                old_axis["time"].min() < log_times.min()
                or old_axis["time"].max() > log_times.max()
            )
        if use_fake_log:
            # TODO ☐: I'm not super interested in considering this type
            #         of file.  If it's the new type, dont' forget that
            #         the log is still there.  If you're talking about
            #         old field/current-sweep experiments, I am fine to
            #         just lok at the DCCT to see how the field changes.
            # These files stored the field directly on the indirect axis, and
            # either have no log or the file-level log is for a different node.
            fake_log = prscr.logobj()
            fake_log.total_log = np.zeros(
                len(old_axis), dtype=fake_log.log_dtype
            )
            fake_log.total_log["time"] = old_axis["time"]
            fake_log.total_log["field"] = old_axis["field"]
            fake_log.total_log["power"] = -np.inf
            fake_log.total_log["cmd"] = 0
            s.set_prop("log", fake_log)
        if vernum < 4:
            # Make fake start/stop times centered on the old experiment time,
            # so we can use the log-averaging code.
            fake_axis = np.zeros(
                len(old_axis),
                dtype=[("start_times", "f8"), ("stop_times", "f8")],
            )
            fake_dt = s.get_prop("acq_params").get("acq_time_ms", 1e3) * 1e-3
            fake_axis["start_times"] = old_axis["time"] - fake_dt / 2
            fake_axis["stop_times"] = old_axis["time"] + fake_dt / 2
            s.setaxis("indirect", fake_axis).set_units("indirect", None)
            s = generate_coordinates_from_log(s, fl=None)
        # }}}

        log_axis = s["indirect"].copy()
        log_axis_error = s.get_error("indirect")

        elapsed = log_axis["time"] - log_axis["time"][0]
        gamma_eff_mhz_g = s.get_prop("acq_params").get(
            "gamma_eff_mhz_g",
            s.get_prop("acq_params").get("gamma_eff_MHz_G"),
        )
        nominal_field_G = (
            s.get_prop("acq_params")["carrierFreq_MHz"] / gamma_eff_mhz_g
        )
        field_shift_Hz = (
            psd.nddata(
                (np.asarray(log_axis["field"]) - nominal_field_G)
                * gamma_eff_mhz_g
                * 1e6,
                [-1],
                ["elapsed"],
            )
            .setaxis("elapsed", elapsed)
            .set_units("elapsed", "s")
        )
        field_shift_error_Hz = log_axis_error["field"] * gamma_eff_mhz_g * 1e6
        s.setaxis("indirect", elapsed).set_units("indirect", "s")
        s.rename("indirect", "elapsed")
        s.ift("t2")
        s *= np.exp(-s.fromaxis("t2") * pi * 50)
        s.ft("t2")
        if "nScans" in s.dimlabels:
            s = s.mean("nScans")
        fl.next("NMR signal - $\\varphi_0$ only")
        s /= prscr.zeroth_order_ph(s)
        fl.DCCT(s)
        fl.next("slice FID")
        s = prscr.fid_from_echo(
            s, signal_pathway=s.get_prop("coherence_pathway")
        )
        s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
        frq_center, frq_half = prscr.find_peakrange(
            s, direct="t2", peak_lower_thresh=0.05
        )
        frq_range = (
            frq_center - frq_half,
            frq_center + frq_half,
        )
        fl.plot(s.real)
        fl.next("NMR signal - with zf and conv (tdom)")
        s.ift("t2").ft("t2", pad=s.shape["t2"] * 20).convolve("t2", 0.25e3)
        fl.plot(s)
        fl.next("B field and peak shift vs. time")
        plt.errorbar(
            elapsed,
            field_shift_Hz.data,
            yerr=field_shift_error_Hz,
            fmt=".",
            label="Hall probe",
        )
        s = s.real["t2":frq_range]
        s /= s.C.integrate("t2")
        fl.next("normalized (zoomed)")
        fl.plot(s)
        s *= s.fromaxis("t2")
        fl.next("freq weighted function")
        fl.plot(s)
        s.integrate("t2")
        fl.next("B field and peak shift vs. time")
        fl.plot(s, "o", label="NMR average freq")
        plt.ylim(frq_range)
        plt.subplots_adjust(hspace=0.5)
        plt.subplots_adjust(wspace=0.3)
