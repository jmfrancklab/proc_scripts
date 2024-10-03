import pyspecdata as psd
import pyspecProcScripts as prscr
import os, h5py
import matplotlib.pyplot as plt
import numpy as np
from Instruments.logobj import logobj
import re

filename, exp_type, nodename, postproc_type = (
    "240924_13p5mM_TEMPOL_ODNP_1.h5",
    "ODNP_NMR_comp/ODNP",
    "ODNP",
    "spincore_ODNP_v5",
)
my_filename = psd.search_filename(
    re.escape(filename), exp_type=expdata, unique=True
)
# {{{ this pulls the log from the file
with h5py.File(my_filename, "r") as f:
    log_grp = f["log"]
    thislog = logobj.from_group(log_grp)
    log_array = thislog.total_log
    log_dict = thislog.log_dict
# }}}
print(log_array.dtype.names)
# I don't know if the following are code is correct, but you should be
# able to correct based on the names of the fields, given by the
# previous line
log_array["time"] -= log_array["time"][0]  # always convert to relative
#                                           time right away
with psd.figlist_var() as fl:
    fl.next("power log")
    log_vs_time = psd.nddata(
        log_array[
            "power"
        ],  # new nddata, whose data are the values from the gigatronix
        [-1],  # it's one dimension, whose length is automatically determined
        ["time"],  # the name of hte dimension is time
    ).setaxis(
        "time",  # set the coordinate axis
        log_array[
            "time"
        ],  # to the "time" field of the structured array that comes from the
        #     log
    )
    fl.plot(log_vs_time)  # should be a picture of the gigatronics powers
    # {{{ construct an nddata whose data are the average power values,
    #     whose errors are the std of of the power values, and whose time
    #     axis is the center time for each power
    # {{{ AG does something else, but basically we want to create an
    #     nddata that will store our powers and the associated errors
    powername = "time"
    dnp_time_axis = s["indirect"]
    log_vs_time.set_units("time", "s")
    mean_power_vs_time = (
        psd.ndshape([("time", len(dnp_time_axis))])
        .alloc()
        .set_error(0)
        .set_units("time", "s")
        .setaxis("time", np.zeros(len(dnp_time_axis)))
    )
    print(dnp_time_axis)
    print("log_vs_time is ", log_vs_time)
    print(log_vs_time.dimlabels)
    relative_times = []
    # }}}
    print(mean_power_vs_time)
    for j, (time_start, time_stop) in enumerate(
        zip(dnp_time_axis[:]["start_times"], dnp_time_axis[:]["stop_times"])
    ):
        print(log_vs_time["time", -1])
        print(log_vs_time["time":(time_start, time_stop)])
        print(mean_power_vs_time[powername, j])
        mean_power_vs_time[powername, j] = log_vs_time[
            "time":(time_start, time_stop)
        ].mean("time", std=True)
        print(mean_power_vs_time[powername, j])
        fl.next("power log")
        relative_times.append(
            (time_stop - log_vs_time.getaxis("time")[0])
            + (time_start - log_vs_time.getaxis("time")[0]) / 2
        )
        # {{{ these lines show the start and the stop of the power step
        plt.axvline(
            x=(time_start - log_vs_time.getaxis("time")[0]),
            color="green",
            alpha=0.5,
        )
        plt.axvline(
            x=(time_stop - log_vs_time.getaxis("time")[0]),
            color="red",
            alpha=0.5,
        )
        # }}}
    print(relative_times)
    mean_power_vs_time.setaxis("time", relative_times)
    fl.plot(
        mean_power_vs_time, "o"
    )  # this  should be a *single* o at the center of each power step.
    #    Its y value should be the avaerage power for that step, and its
    #    error bars should give the standard deviation of the power over
    #    the step
    # }}}
fl.show()
