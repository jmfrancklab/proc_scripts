import pyspecdata as psd
import h5py
import matplotlib.pyplot as plt
import numpy as np
from Instruments.logobj import logobj
import datetime
import re


def plot_field(
    s,
    filename,
    exp_type,
    fl=None,
    node_name="log",
):
    """Generate field axis for ODNP/E(p)/FIR experiments, converts instrument
    field log to useable axis. Adapted from convert_to_power.py

    Parameters
    ===========
    s: nddata
        This is an nddata whose "indirect"
        axis is a structured array.  This structured array gives the
        `start_times` and `stop_times` (as fields) for each datapoint
        along the indirect axis.
    filename : str
        Name of the hdf file that contains the field log.
    exp_type : str
        The "experiment type" used to search for the file.
    fl: figlist_var()

    Returns
    =======
    values : numpy.array
    s : nddata
        (modified in-place)
        When we're done, the axis coordinates for "indirect" are a (not
        structured) array containing the average field recorded on the
        log between each start_time and stop_time.
        The axis coordinate errors are the standard deviations of the
        same.
    """
    my_filename = psd.search_filename(
        re.escape(filename), exp_type=exp_type, unique=True
    )
    with h5py.File(my_filename, "r") as f:
        thislog = logobj.from_group(f[node_name])
        log_array = np.array(thislog.total_log)

    assert "field" in log_array.dtype.names, str(log_array.dtype.names)
    assert s["indirect"].dtype.names == (
        "start_times",
        "stop_times",
    ), str(s["indirect"].dtype.names)

    zero_time = log_array["time"][0]
    log_array["time"] -= zero_time

    log_vs_time = (
        psd.nddata(
            log_array["field"],  # Create an nddata object from the field data
            [
                -1
            ],  # The shape of the data is 1D, so we use -1 to indicate that the length should be inferred
            ["time"],  # The axis is named "time"
        )
        .setaxis("time", log_array["time"])
        .set_units("time", "s")
    )
    log_vs_time.set_units("G")

    if fl:  # Check that fl is not None before plotting
        fl.next("field log")
        fl.plot(
            log_vs_time,
            ".",
            human_units=False,
        )  # Plot the field log data as points
        ax = plt.gca()
        # {{{ this is just matplotlib time formatting
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(
                lambda x, _: (
                    str(datetime.timedelta(seconds=x)).lstrip("0:").lstrip(":")
                    if x > 0
                    else "0:00"
                )
            )
        )
        # }}}
    log_vs_time.set_units("time", "s")
    mean_field_vs_time = (
        psd.ndshape([("time", len(s["indirect"]))])
        .alloc(dtype=np.float64)
        .set_error(0)
        .set_units("time", "s")
        .setaxis("time", np.zeros(len(s["indirect"])))
    )
    # {{{ we need to convert these to relative time up front, so that things
    #     don't get complicated!
    s["indirect"]["start_times"] -= zero_time
    s["indirect"]["stop_times"] -= zero_time
    # }}}
    for j, (time_start, time_stop) in enumerate(
        zip(
            s["indirect"][:]["start_times"],
            s["indirect"][:]["stop_times"],
        )
    ):
        mean_field_vs_time["time", j] = log_vs_time[
            "time" : (time_start, time_stop)
        ].mean("time", std=True)
        mean_field_vs_time["time"][j] = (time_start + time_stop) / 2
        plt.axvspan(
            time_start,
            time_stop,
            facecolor="none",
            edgecolor="k",
            hatch="XXXXXX",
            alpha=0.1,
        )
        mean_field_vs_time.set_units("G")
    # {{{ Plot mean field vs. time with error bars
    if fl:
        fl.plot(
            mean_field_vs_time,
            "o",
            human_units=False,
        )

    s.setaxis("indirect", mean_field_vs_time.data).set_error(
        "indirect", mean_field_vs_time.get_error()
    ).set_units("indirect", "G")
    s.rename("indirect", "field")
    # }}}
    return s
