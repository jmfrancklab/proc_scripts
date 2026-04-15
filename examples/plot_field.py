# TODO ☐: not yet reviewed

import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt


# {{{ the HDF used in this particular example is broken, so we need to
#     patch it
def _decode_list_node(h5group):
    item_names = sorted(
        (name for name in h5group.attrs if name.startswith("ITEM")),
        key=lambda name: int(name[4:]),
    )
    values = []
    for name in item_names:
        value = h5group.attrs[name]
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        values.append(value)
    return values


def fix_broken_hdf(log_group):
    # {{{ because this is a hack, let's just create our classes inline,
    #     to keep it simple
    array_node_cls = type(
        "BrokenArrayNode",
        (),
        {
            "__getitem__": lambda self, item: self._array[item],
        },
    )
    group_node_cls = type(
        "BrokenGroupNode",
        (),
        {
            "keys": lambda self: ["array"],
            "__getitem__": lambda self, key: (
                self._array_node
                if key == "array"
                else (_ for _ in ()).throw(KeyError(key))
            ),
        },
    )
    # }}}
    array_node = array_node_cls()
    array_node._array = log_group["array"][:]
    array_node.attrs = {
        "dictkeys": _decode_list_node(log_group["dictkeys"]),
        "dictvalues": _decode_list_node(log_group["dictvalues"]),
    }
    group_node = group_node_cls()
    group_node._array_node = array_node
    return group_node


# }}}


filename = "260406_hydroxytempo_ODNP_1.h5"
exp_type = "B27/ODNP"
# TODO: I'm not stoked about this -- can't we just search_filename, and
#       then use the load_indiv? or, alternately (would be better for
#       this log stuff), should we modify
#       find_file to return the full local path as a property?
#       Even better, should we modify the functions in the lookup table
#       so that, for the right postproc types, it goes and grabs the log
#       and attaches it as a property.
log_array = prscr.load_log_data(filename, exp_type, hdf_repair=fix_broken_hdf)
s = psd.find_file(
    filename,
    exp_type=exp_type,
    expno="ODNP",
    lookup=prscr.lookup_table,
)
zero_time = log_array["time"][0].item()
log_array["time"] -= zero_time
carrier_Hz = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
gamma_eff_Hz_G = s.get_prop("acq_params")["gamma_eff_MHz_G"] * 1e6
field_drift_Hz = log_array["field"] * gamma_eff_Hz_G - carrier_Hz
s.rename("indirect", "t").set_units("t", "s")
s["t"] = 0.5 * (s["t"]["start_times"] + s["t"]["stop_times"]) - zero_time
# {{{ rather than averaging over nScans, I create an appropriate time
#     axis for it, and then smoosh to create a dimension with all my scans
s["nScans"] = s["nScans"] * (
    s.get_prop("acq_params")["acq_time_ms"] * 1e-3
    + s.get_prop("acq_params")["repetition_us"] * 1e-6
)
s.smoosh(["t", "nScans"], dimname="t")
# smoosh makes a record array
s["t"] = s["t"]["nScans"] + s["t"]["t"]
s.set_units("t", "s")
# }}}
s.reorder("t2", first=False)
field_drift_Hz = (
    psd.nddata(field_drift_Hz, [-1], ["t"])
    .setaxis("t", log_array["time"])
    .set_units("t", "s")
)
with psd.figlist_var() as fl:
    fl.next("NMR signal - $\\varphi_0$ only")
    s /= prscr.zeroth_order_ph(s)
    fl.DCCT(s)
    fl.next("slice FID")
    s = prscr.fid_from_echo(s, signal_pathway=s.get_prop("coherence_pathway"))
    s = prscr.select_pathway(s, s.get_prop("coherence_pathway"))
    frq_range = (-5e3, 5e3)
    fl.plot(s.real)
    fl.next("NMR signal - with zf and conv (tdom)")
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
    fl.plot(field_drift_Hz, ".", label="Hall probe")
    fl.plot(s, "o", label="NMR")
    plt.ylim(frq_range)
