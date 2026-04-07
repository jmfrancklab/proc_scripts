import pyspecdata as psd
import pyspecProcScripts as prscr

filename = "260406_hydroxytempo_ODNP_1.h5"
exp_type = "B27/ODNP"
# TODO: I'm not stoked about this -- can't we just search_filename, and
#       then use the load_indiv? or, alternately (would be better for
#       this log stuff), should we modify
#       find_file to return the full local path as a property?
#       Even better, should we modify the functions in the lookup table
#       so that, for the right postproc types, it goes and grabs the log
#       and attaches it as a property.
log_array = prscr.load_log_data(filename, exp_type)
s = psd.find_file(
    filename,
    exp_type=exp_type,
    expno="ODNP",
    lookup=prscr.lookup_table,
)
s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
if "nScans" in s.dimlabels:
    s = s.mean("nScans")
zero_time = log_array["time"][0]
log_array["time"] -= zero_time
s["indirect"] -= zero_time
carrier_Hz = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
gamma_eff_Hz_G = s.get_prop("acq_params")["gamma_eff_MHz_G"] * 1e6
field_drift_Hz = log_array["field"] * gamma_eff_Hz_G - carrier_Hz
s.rename("indirect", "t").set_units("t", "s")
s["t"] = 0.5 * (s["indirect"]["start_times"] + s["indirect"]["stop_times"])
field_drift_Hz = (
    psd.nddata(field_drift_Hz, [-1], ["t"])
    .setaxis("t", log_array["t"])
    .set_units("t", "s")
)
with psd.figlist_var() as fl:
    fl.next("raw NMR signal")
    fl.image(s)
    s.run(abs).argmax("t2")
    fl.next("B field and peak shift vs. time")
    fl.plot(field_drift_Hz, ".", label="Hall probe")
    fl.plot(s, "o", label="NMR")
