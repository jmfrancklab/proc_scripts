import datetime

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

import pyspecProcScripts as prscr
from .generate_coordinates_from_log import load_log_data


def attach_log_data_from_file(
    s,
    filename,
    exp_type,
    node_name="log",
):
    """Attach legacy HDF log data as an nddata property."""
    log_array = load_log_data(
        filename,
        exp_type,
        node_name=node_name,
    )
    thislog = prscr.logobj()
    thislog.total_log = log_array
    s.set_prop("log", thislog)
    return s


def generate_power_coordinates_from_log(
    s,
    filename,
    exp_type,
    fl=None,
    node_name="log",
    directional_coupler_dB=22,
):
    """Generate a scalar power axis from HDF log data.

    This preserves the old convert_to_power behavior for pre-v5 postproc
    files, where the log is not attached to the nddata by the loader.
    """
    s = attach_log_data_from_file(
        s,
        filename,
        exp_type,
        node_name=node_name,
    )
    log_array = s.get_prop("log").total_log.copy()
    assert all(
        name in log_array.dtype.names
        for name in ["time", "Rx", "power", "cmd"]
    ), str(log_array.dtype.names)
    assert s["indirect"].dtype.names == (
        "start_times",
        "stop_times",
    ), str(s["indirect"].dtype.names)
    zero_time = log_array["time"][0]
    log_array["time"] -= zero_time
    log_vs_time = (
        psd.nddata(log_array["power"], [-1], ["time"])
        .setaxis("time", log_array["time"])
        .set_units("time", "s")
    )
    log_vs_time = prscr.dBm2power(log_vs_time + directional_coupler_dB)
    log_vs_time.set_units("W")
    if fl:
        fl.next("power log")
        fl.plot(log_vs_time, ".", human_units=False)
        ax = plt.gca()
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(
                lambda x, _: (
                    str(datetime.timedelta(seconds=x)).lstrip("0:").lstrip(":")
                    if x > 0
                    else "0:00"
                )
            )
        )
    mean_power_vs_time = (
        psd.ndshape([("time", len(s["indirect"]))])
        .alloc(dtype=np.float64)
        .set_error(0)
        .set_units("time", "s")
        .setaxis("time", np.zeros(len(s["indirect"])))
    )
    s["indirect"]["start_times"] -= zero_time
    s["indirect"]["stop_times"] -= zero_time
    for j, (time_start, time_stop) in enumerate(
        zip(
            s["indirect"][:]["start_times"],
            s["indirect"][:]["stop_times"],
        )
    ):
        mean_power_vs_time["time", j] = log_vs_time[
            "time":(time_start, time_stop)
        ].mean("time", std=True)
        mean_power_vs_time["time"][j] = (time_start + time_stop) / 2
        plt.axvspan(
            time_start,
            time_stop,
            facecolor="none",
            edgecolor="k",
            hatch="XXXXXX",
            alpha=0.1,
        )
        mean_power_vs_time.set_units("W")
    if fl:
        fl.plot(mean_power_vs_time, "o", human_units=False)
    s.setaxis("indirect", mean_power_vs_time.data).set_error(
        "indirect", mean_power_vs_time.get_error()
    ).set_units("indirect", "W")
    s["indirect"][abs(s["indirect"]) < 10**-10] = 0
    s.rename("indirect", "power")
    return s
