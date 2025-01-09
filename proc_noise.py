from pylab import *
from pyspecdata import *
from scipy.optimize import leastsq, minimize, basinhopping
from itertools import cycle
import matplotlib.pyplot as plt

fl = figlist_var()
scalefactor = 35000  # a manual adjustment -- seems to bring the signal level to close to 1
convwidth = 0.5e3
frq_slice = (-50e3, 50e3)
showraw = False
new_colors = cycle(
    [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
)
for date, id_string, label_str in [
    ("200210", "echo_SW_11", "SW=3 kHz"),
    ("200210", "echo_SW_10", "SW=6 kHz"),
    ("200210", "echo_SW_9", "SW=12 kHz"),
    ("200210", "echo_SW_1", "SW=24 kHz"),
    ("200210", "echo_SW_2", "SW=48 kHz"),
    ("200210", "echo_SW_3", "SW=96 kHz"),
    ("200210", "echo_SW_4", "SW=192 kHz"),
    ("200210", "echo_SW_5", "SW=384 kHz"),
    ("200210", "echo_SW_6", "SW=768 kHz"),
]:
    filename = date + "_" + id_string + ".h5"
    nodename = "signal"
    s = nddata_hdf5(
        filename + "/" + nodename,
        directory=getDATADIR(exp_type="ODNP_NMR_comp/old"),
    )
    s.name(None)
    s /= scalefactor
    nPoints = s.get_prop("acq_params")["nPoints"]
    nEchoes = s.get_prop("acq_params")["nEchoes"]
    nPhaseSteps = s.get_prop("acq_params")["nPhaseSteps"]
    SW_kHz = s.get_prop("acq_params")["SW_kHz"]
    nScans = s.get_prop("acq_params")["nScans"]
    assert nPhaseSteps == 8
    s.chunk("t", ["ph2", "ph1", "t2"], [2, 4, -1])
    s.setaxis("ph2", r_[0.0, 2.0] / 4)
    s.setaxis("ph1", r_[0.0, 1.0, 2.0, 3.0] / 4)
    s.setaxis("nScans", r_[0:nScans])
    s.set_units("t2", "s")
    signal_len = s.getaxis("t2")[-1]
    s.ft("t2", shift=True)
    s.reorder(["ph1", "ph2", "nScans"])
    if showraw:
        fl.basename = label_str
        fl.next("raw data, chunked")
        fl.image(s)
    noise = s.C.run(abs).run(lambda x: x**2).mean_all_but("t2")
    noise /= signal_len  # convert to POWER
    s.ft(["ph1", "ph2"])
    if showraw:
        fl.next("coherence")
        fl.image(s)
        fl.basename = None
    fl.next("signal", legend=True)
    signal = s["t2":(-1000, 1000)]["ph1", 1]["ph2", -2]
    c = next(new_colors)
    signal.ift("t2")
    fl.plot(
        abs(signal)["t2":(None, 0.2)],
        linewidth=3,
        color=c,
        alpha=0.25,
        label=label_str,
        human_units=False,
    )  # here, it makes sense to set
    #                     human_units to false because the spectra are very
    #                     different widths -- because human_units calculates
    #                     the default units each line independently, and sees
    #                     if the new line falls in with previous, this will not
    #                     work
    fl.next("noise -- semilog", legend=True)
    fl.plot(
        noise["t2":frq_slice],
        plottype="semilogy",
        alpha=0.05,
        color=c,
        human_units=False,
    )
    fl.plot(
        noise["t2":frq_slice],
        plottype="semilogy",
        alpha=0.5,
        color=c,
        label=label_str,
        human_units=False,
    )
    fl.next("noise -- linear", legend=True)
    fl.plot(
        noise["t2":frq_slice],
        alpha=0.5,
        color=c,
        label=label_str,
        human_units=False,
    )
    plt.ylim(0, 1e-4)
fl.show()
