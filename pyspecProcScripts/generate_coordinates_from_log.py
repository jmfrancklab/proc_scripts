import pyspecdata as psd
import pyspecProcScripts as prscr
import matplotlib.pyplot as plt
import numpy as np
import datetime


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
        When we're done, the axis coordinates for "indirect" are a
        structured array containing the original time coordinate together
        with the averages of the logged quantities recorded on the
        log between each start_time and stop_time.
        The axis coordinate errors are the standard deviations of the
        same.
    """
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
    log_array["power"] = prscr.dBm2power(
        log_array["power"] + directional_coupler_dB
    )
    log_vs_time = (
        psd.nddata(log_array, [-1], ["time"])
        .set_axis("time", log_array["time"])
        .set_units("time", "s")
    )
    if fl:  # checks that fl is not None
        fl.next("power log")
        # TODO ☐: I know that this doesn't work.  Use a loop to plot
        #         over subplots.  **Reuse** code that you already wrote
        #         for the read log example that does this.
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
    # TODO ☐: the name of this is bad -- it should be the mean of
    #         everything, right?
    mean_log_quant_vs_time = (
        psd.ndshape([("time", len(s["indirect"]))])
        .alloc(dtype=np.float64)
        .set_error(0)
        .set_units("time", "s")
        .set_axis("time", np.zeros(len(s["indirect"])))
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
        mean_log_quant_vs_time["time", j] = log_vs_time[
            "time" : (time_start, time_stop)
        ].mean("time", std=True)
        # {{{ I realized a crosshatch would be better here
        plt.axvspan(
            time_start,
            time_stop,
            facecolor="none",
            edgecolor="k",
            hatch="XXXXXX",
            alpha=0.1,
        )
        # mean_log_quant_vs_time = prscr.dBm2power(mean_log_quant_vs_time)
        # }}}
    if fl:
        fl.plot(
            mean_log_quant_vs_time,
            "o",
            human_units=False,
        )  # this  should be a *single* o at the center of each power step.
        #    Its y value should be the avaerage power for that step, and its
        #    error bars should give the standard deviation of the power over
        #    the step
    # }}}
    s.setaxis("indirect", mean_log_quant_vs_time.data).set_error(
        "indirect", mean_log_quant_vs_time.get_error()
    ).set_units("indirect", None) # for now, we need to set this to no units
    s["indirect"][abs(s["indirect"]) < 10**-10] = 0  # the power log
    #                                                 reads as a very
    #                                                 very small power
    #                                                 rather than 0, so
    #                                                 threshold these
    #                                                 out
    return s
