import pyspecdata as psd
import pyspecProcScripts as prscr
import os, h5py
import matplotlib.pyplot as plt
import numpy as np
from Instruments.logobj import logobj
import matplotlib.dates as mdates
import datetime
import re

filename, exp_type, nodename, postproc_type = (
    "240924_13p5mM_TEMPOL_ODNP_1.h5",
    "ODNP_NMR_comp/ODNP",
    "ODNP",
    "spincore_ODNP_v5",
)
# {{{ we need to manually pull the log (a structured array and a dictionary)
#     and to pull the spectrum `s`, which is an nddata
my_filename = psd.search_filename(
    re.escape(filename), exp_type=exp_type, unique=True
)
with h5py.File(my_filename, "r") as f:
    thislog = logobj.from_group(f["log"])
    log_array = thislog.total_log
    log_dict = thislog.log_dict
assert log_array.dtype.names == ("time", "Rx", "power", "cmd"), str(
    log_array.dtype.names
)
# I don't know if the following are code is correct, but you should be
# able to correct based on the names of the fields, given by the
# previous line
zero_time = log_array["time"][0]
log_array["time"] -= zero_time  # always convert to relative time
#                                            right away
s = psd.find_file(re.escape(filename), exp_type=exp_type, expno=nodename, lookup=prscr.lookup_table)
# }}}
with psd.figlist_var() as fl:
    fl.next("power log")
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
            ],  # to the "time" field of the structured array that comes from the
            #     log
        )
        .set_units("time", "s")
    )
    fl.plot(log_vs_time, ".")  # should be a picture of the gigatronics powers
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
    powername = "time"
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
        zip(s["indirect"][:]["start_times"], s["indirect"][:]["stop_times"])
    ):
        mean_power_vs_time[powername, j] = log_vs_time[
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
        # }}}
    fl.plot(
        mean_power_vs_time, "o"
    )  # this  should be a *single* o at the center of each power step.
    #    Its y value should be the avaerage power for that step, and its
    #    error bars should give the standard deviation of the power over
    #    the step
    # }}}
