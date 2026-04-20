"""
Read Instrument Log
===================
"""

from pyspecProcScripts import logobj
import pyspecdata as psd
import h5py
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import datetime


def fix_broken_hdf(log_group):
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


coupler_atten = 22
files_to_check = [
    ("230626_batch230515_E37_Ras_B10_ODNP_1.h5", "ODNP_NMR_comp/ODNP", False),
    ("260406_hydroxytempo.*", "B27/ODNP", True),
    ("260107_hydroxytempo_ODNP_1.h5", "B27/ODNP", False),
]

with psd.figlist_var() as fl:
    for which_fname, which_exptype, broken_hdf in files_to_check:
        myfilename = psd.search_filename(
            which_fname,
            exp_type=which_exptype,
            unique=True,
        )
        # {{{ open h5 file to real log
        with h5py.File(myfilename, "r") as f:
            if broken_hdf:
                thislog = logobj.from_group(fix_broken_hdf(f["log"]))
            else:
                thislog = logobj.from_group(f["log"])
        # }}}
        # In order to properly set the time axis to start at 0
        # both the log's start time will be subtracted from the
        # the relative time recorded
        print("ask for overall type", type(thislog.total_log))
        print("ask for dtype", thislog.total_log.dtype)
        thislog.total_log["time"] -= thislog.total_log["time"][0]
        plot_field = "field" in thislog.total_log.dtype.names
        # }}}
        # {{{ plot the output power and reflection
        fig, ax_list = plt.subplots(3 if plot_field else 2, 1, figsize=(10, 8))
        fl.next("log figure", fig=fig)
        if plot_field:
            ax_Rx, ax_power, ax_field = ax_list
            ax_field.set_ylabel("field / G")
            ax_field.set_xlabel("Time / ms")
            ax_field.plot(
                thislog.total_log["time"],
                thislog.total_log["field"],
                "k,",
                alpha=0.25,
            )
        else:
            ax_Rx, ax_power = ax_list
        ax_Rx.set_ylabel("Rx / mV")
        ax_Rx.set_xlabel("Time / ms")
        ax_Rx.plot(
            thislog.total_log["time"], thislog.total_log["Rx"], "k,", alpha=0.25
        )
        ax_power.set_ylabel("power / dBm")
        ax_power.set_xlabel("Time / ms")
        ax_power.plot(
            thislog.total_log["time"],
            10
            ** (
                (thislog.total_log["power"] + coupler_atten) / 10 - 3
            ),  # -3 for mW to W
            "k,",
            alpha=0.25,
        )
        # }}}
        mask = thislog.total_log["cmd"] != 0
        position = 0
        npositions = 20
        for j, thisevent in enumerate(thislog.total_log[mask]):
            # {{{ Add a vertical line at the time the data acquisition for the
            #     set power began
            print(
                "looking for key",
                thisevent["cmd"],
                "in dictionary with keys",
                thislog.log_dict.keys(),
            )
            event_name = thislog.log_dict[thisevent["cmd"]]
            if event_name.lower().startswith("get_power"):
                continue  # ignore "get power" commands
            position = (
                position % npositions
            )  # use npositions positions top to bottom, then roll over
            for thisax in ax_list:
                thisax.axvline(x=thisevent["time"], color="g", alpha=0.3)
                thisax.text(
                    s=event_name,
                    x=thisevent["time"],
                    y=0.1 + (0.9 - 0.1) * position / npositions,
                    transform=blended_transform_factory(
                        thisax.transData, thisax.transAxes
                    ),
                    alpha=0.3,
                    color="g",
                    size=8,  # really tiny!
                )
            position += 1
            # }}}
        for thisax in ax_list:
            thisax.xaxis.set_major_formatter(
                plt.FuncFormatter(
                    lambda x, _: str(datetime.timedelta(seconds=x))
                )
            )
        plt.tight_layout()
