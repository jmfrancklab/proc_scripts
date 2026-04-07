import pyspecdata as psd
import pyspecProcScripts as prscr

filename = "260406_hydroxytempo_ODNP_1.h5"
exp_type = "B27/ODNP"
nodename = "ODNP"

log_array = prscr.load_log_data(filename, exp_type)
s = psd.find_file(
    filename, exp_type=exp_type, expno=nodename, lookup=prscr.lookup_table
)
s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
if "nScans" in s.dimlabels:
    s = s.mean("nScans")
if not s.get_ft_prop("t2"):
    s.ft("t2")
zero_time = log_array["time"][0]
log_array["time"] -= zero_time
carrier = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
gamma_eff_Hz_G = s.get_prop("acq_params")["gamma_eff_MHz_G"] * 1e6
field_shift = log_array["field"] * gamma_eff_Hz_G - carrier
time_midpoints = (
    0.5 * (s["indirect"]["start_times"] + s["indirect"]["stop_times"])
    - zero_time
)
field_shift_vs_time = (
    psd.nddata(field_shift, [-1], ["time"])
    .setaxis("time", log_array["time"])
    .set_units("time", "s")
)
peak_shift_vs_time = s.run(abs).argmax("t2")
peak_shift_vs_time.rename("indirect", "time")
peak_shift_vs_time.setaxis("time", time_midpoints).set_units("time", "s")

with psd.figlist_var() as fl:
    fl.next("B field and peak shift vs. time")
    fl.plot(field_shift_vs_time, ".", label="field shift")
    fl.plot(peak_shift_vs_time, "o", label="peak shift")
