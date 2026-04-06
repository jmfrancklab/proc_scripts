"""
Plot raw logged power against time
=================================
Read the instrument log directly and plot the recorded power as a function of
time without using convert_to_power().
"""

import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt

thisfile, exp_type, nodename = (
    "260107_hydroxytempo_ODNP_1.h5",
    "B27/ODNP",
    "ODNP",
)
log_array = prscr.load_log_data(thisfile, exp_type)
log_array["time"] -= log_array["time"][0]
power_vs_time = (
    psd.nddata(log_array["power"], [-1], ["time"])
    .setaxis("time", log_array["time"])
    .set_units("time", "s")
)
data = psd.find_file(
    thisfile,
    exp_type=exp_type,
    expno=nodename,
    lookup=prscr.lookup_table,
)
data = prscr.select_pathway(data, data.get_prop("coherence_pathway"))
orig_axis = data["indirect"]
print(orig_axis["start_times"], orig_axis["stop_times"])
acq_times = -orig_axis["start_times"] + orig_axis["stop_times"]

print(acq_times)
acq_times -= acq_times[0, :]
data["indirect"] = acq_times
data.set_units("indirect", "s")
forplot = data.run(abs).argmax("t2")
with psd.figlist_var() as fl:
    fl.next("power and peak vs time")
    ax_peak = plt.gca()
    fl.plot(forplot, "o", label="peak", color="C0")
    ax_peak.set_ylabel("peak", color="C0")
    ax_peak.tick_params(axis="y", colors="C0")
    ax_peak.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0))
    ax_power = ax_peak.twinx()
    fl.plot(power_vs_time, ".", ax=ax_power, label="power", color="C1")
    ax_power.set_ylabel("power", color="C1")
    ax_power.tick_params(axis="y", colors="C1")
    ax_power.legend(loc="lower right", bbox_to_anchor=(1.0, 0.9))
