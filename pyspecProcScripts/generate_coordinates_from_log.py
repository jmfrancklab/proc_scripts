import pyspecdata as psd
import pyspecProcScripts as prscr
import h5py
import matplotlib.pyplot as plt
import numpy as np
import datetime
import re


def load_log_data(
    filename,
    exp_type,
    node_name="log",
    hdf_repair=None,
):
    """Load instrument log from an HDF5 file.

    Parameters
    ==========
    hdf_repair: function default None
        For some intermediate versions with broken HDF storage, this
        allows us to supply a patch function that fixes the data.
    """
    filename = psd.search_filename(
        re.escape(filename), exp_type=exp_type, unique=True
    )
    with h5py.File(filename, "r") as f:
        if hdf_repair is None:
            thislog = prscr.logobj.from_group(f[node_name])
        else:
            thislog = prscr.logobj.from_group(hdf_repair(f[node_name]))
        log_array = np.array(thislog.total_log, copy=True)
    return log_array


def generate_coordinates_from_log(
    s,
    filename,
    exp_type,
    fl=None,
    node_name="log",
    directional_coupler_dB=22,
):
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
        directional_coupler_dB : float
        The ratio (in dB) of the power that leaves the main port of the
        directional coupler vs. the power that arrives at the power
        detector.  This is measured, and the value here is a previously
        measured value for our system.
    fl: figlist_var()

    Returns
    =======

    values : numpy.array
    s : nddata
        (modified in-place)
        When we're done, the axis coordinates for "indirect" are a (not
        structured) array containing the average power recorded on the
        log between each start_time and stop_time.
        The axis coordinate errors are the standard deviations of the
        same.
    """
    log_array = s.get_prop("log")
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
    log_vs_time["power"] = prscr.dBm2power(
        log_vs_time["power"] + directional_coupler_dB
    )
    if fl:  # checks that fl is not None
        fl.next("power log")
        fl.plot(
            log_vs_time,
            ".",
            human_units=False,
        )  # should be a picture of the gigatronics powers
        # {{{ this is just matplotlib time formatting
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
        # }}}
    # {{{ construct an nddata whose data are the average power values,
    #     whose errors are the std of of the power values, and whose time
    #     axis is the center time for each power
    # {{{ AG does something else, but basically we want to create an
    #     nddata that will store our powers and the associated errors
    log_vs_time.set_units("time", "s")
    mean_power_vs_time = (
        psd.ndshape([("time", len(s["indirect"]))])
        .alloc(dtype=np.float64)
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
        mean_power_vs_time["time", j] = log_vs_time.mean("time", std=True)
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

        # }}}
    if fl:
        fl.plot(
            mean_power_vs_time,
            "o",
            human_units=False,
        )  # this  should be a *single* o at the center of each power step.
        #    Its y value should be the avaerage power for that step, and its
        #    error bars should give the standard deviation of the power over
        #    the step
    # }}}
    s.setaxis("indirect", mean_power_vs_time.data).set_error(
        "indirect", mean_power_vs_time.get_error()
    )
    s["indirect"]["power"][abs(s["indirect"]) < 10**-10] = 0  # the power log
    #                                                 reads as a very
    #                                                 very small power
    #                                                 rather than 0, so
    #                                                 threshold these
    #                                                 out
    return s
