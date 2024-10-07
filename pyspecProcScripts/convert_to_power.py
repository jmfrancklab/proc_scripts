import pyspecdata as psd
import pyspecProcScripts as prscr
import h5py
import matplotlib.pyplot as plt
import numpy as np
from Instruments.logobj import logobj
import datetime
import re


def convert_to_power(s, filename, exp_type, fl=None, node_name="log"):
    """Generate power axis for ODNP/E(p)/FIR experiments, converts instrument
    power log to useable axis
    Parameters
    ===========
    s: nddata
        This is an nddata whose "indirect"
        axis is a structured array.  This structured array gives the
        `start_times` and `stop_times` (as fields) for each datapoint
        along the indirect axis.
    filename : str
        Name of the hdf file that contains the power log.
    exp_type : str
        The "experiment type" used to search for the file.
    fl: figlist_var()

    Returns
    =======
    values: numpy.array
    s: nddata
        (modified in-place)
        When we're done, the axis coordinates for "indirect" are a (not
        structured) array containing the average power recorded on the
        log between each start_time and stop_time.
        The axis coordinate errors are the standard deviations of the
        same.
    """
    my_filename = psd.search_filename(
        re.escape(filename), exp_type=exp_type, unique=True
    )
    with h5py.File(my_filename, "r") as f:
        thislog = logobj.from_group(f[node_name])
        log_array = thislog.total_log
    assert log_array.dtype.names == ("time", "Rx", "power", "cmd"), str(
        log_array.dtype.names
    )
    assert s["indirect"].dtype.names == (
        "start_times",
        "stop_times",
    ), str(s["indirect"].dtype.names)
    zero_time = log_array["time"][0]
    log_array["time"] -= zero_time
    log_vs_time = (
        psd.nddata(
            log_array[
                "power"
            ],  # new nddata, whose data are the values from the gigatronix
            [
                -1
            ],  # it's one dimension, whose length is automatically determined
            ["time"],  # the name of hte dimension is time
        )
        .setaxis(
            "time",  # set the coordinate axis
            log_array[
                "time"
            ],  # to the "time" field of the structured array that comes from
            #     the log
        )
        .set_units("time", "s")
    )
    # TODO ☑: these are the power readings on the gigatronics.  They
    # need to be converted to actual powers by adding in a constant dB
    # value that represents our directional coupler and the associated
    # losses.  This number is about (but not exactly) 20.  Alex should
    # be able to tell you what it is.
    log_vs_time = prscr.dBm2power(log_vs_time + 22)
    log_vs_time.set_units("W")
    if fl:  # checks that fl is not None
        fl.next("power log")
        fl.plot(
            log_vs_time, "."
        )  # should be a picture of the gigatronics powers
        # {{{ this is just matplotlib time formatting
        ax = plt.gca()
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(
                lambda x, _: str(datetime.timedelta(seconds=x))
                .lstrip("0:")
                .lstrip(":")
                if x > 0
                else "0:00"
            )
        )
        # }}}
    # {{{ construct an nddata whose data are the average power values,
    #     whose errors are the std of of the power values, and whose time
    #     axis is the center time for each power
    # {{{ AG does something else, but basically we want to create an
    #     nddata that will store our powers and the associated errors
    log_vs_time.set_units("time", "s")
    mean_power_vs_time = (
        psd.ndshape([("time", len(s["indirect"]))])
        .alloc()
        .set_error(0)
        .set_units("time", "s")
        .setaxis("time", np.zeros(len(s["indirect"])))
    )
    # }}}
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
        mean_power_vs_time["time", j] = log_vs_time[
            "time":(time_start, time_stop)
        ].mean("time", std=True)
        mean_power_vs_time["time"][j] = (time_start + time_stop) / 2
        # {{{ I realized a crosshatch would be better here
        plt.axvspan(
            time_start,
            time_stop,
            facecolor="none",
            edgecolor="k",
            hatch="XXXXXX",
            alpha=0.1,
        )
        # mean_power_vs_time = prscr.dBm2power(mean_power_vs_time)
        mean_power_vs_time.set_units("W")
        # }}}
    if fl:
        fl.plot(
            mean_power_vs_time, "o"
        )  # this  should be a *single* o at the center of each power step.
        #    Its y value should be the avaerage power for that step, and its
        #    error bars should give the standard deviation of the power over
        #    the step
    # }}}
    s.setaxis("indirect", mean_power_vs_time.data).set_error(
        "indirect", mean_power_vs_time.get_error()
    ).set_units("indirect", "W")
    return s
