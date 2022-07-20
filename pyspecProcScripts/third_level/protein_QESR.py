import numpy as np
from numpy import sqrt, r_
import os
from matplotlib.pyplot import axvline, axhline, gca
import matplotlib.pyplot as plt
from pint import UnitRegistry
from pyspecdata.datadir import pyspec_config # piggyback on _pyspecdata
from pyspecdata import find_file, gammabar_e, strm, ndshape
from scipy.interpolate import UnivariateSpline
import pickle

# {{{ all basic info
default_Q = float(pyspec_config.get_setting("default Q", default=4700))
# note that this will fail if you don't have "QESR Conversion" set in your pyspecdata file
dint_conversion = float(pyspec_config.get_setting("QESR Conversion"))

ureg = UnitRegistry(
    system="mks", autoconvert_offset_to_baseunit=True, auto_reduce_dimensions=True
)
Q_ = ureg.Quantity
colors = plt.rcParams[
    "axes.prop_cycle"
]()  # this is the default matplotlib cycler for line styles
fieldaxis = "$B_0$"
#pushout = 1
myconcs = []
# }}}

def protein_QESR(file_name, label, pushout=0.5, threshold=0.05, pickle_file=None, background=None, fl=None, exp_type="francklab_esr/Farhana", which_plot=None):
    if which_plot is None:
        which_plot = file_name
    if fl is None:
        raise ValueError("for now, you just have to pass a figure list")
    if pickle_file is not None:
        # {{{ make a pickle file, all_concs.pickle that has a
        #     dictionary of the concentrations of all samples
        #     that we are looking at
        #     the label should have the batch date and
        #     residue, so should be good identifier
        if os.path.exists(pickle_file):
            with open(pickle_file, "rb") as fp:
                vars = pickle.load(fp)
        else:
            vars = {}
        # }}}
    # {{{ load the file of interest, and get set up
    d = find_file(file_name, exp_type=exp_type)
    if "harmonic" in d.dimlabels:
        d = d["harmonic", 0]
    d -= d[fieldaxis, -100:].data.mean()
    d.setaxis(fieldaxis,
            lambda x: x-d.get_prop('MWFQ')/gammabar_e*1e4)
    # }}}
    if background is None:
        background = ndshape(d).alloc() # zeros -- keep life easy!
        background.copy_axes(d)
    else:
        # {{{ set up the background
        background = background.C
        background -= background[fieldaxis, -100:].data.mean()
        background.setaxis(fieldaxis,
                lambda x: x-background.get_prop('MWFQ')/gammabar_e*1e4)
        # }}}
    background = background.interp(fieldaxis, d.getaxis(fieldaxis))
    # {{{ determine the signal denominator from the parameters of interest
    G_R = Q_(*d.get_prop("Gain"))
    C_t = Q_(*d.get_prop("ConvTime"))
    power = Q_(*d.get_prop("Power"))
    B_m = Q_(*d.get_prop("ModAmp"))
    Q = Q_(default_Q, "dimensionless")  # hard set Q value
    n_B = Q_(1, "dimensionless")  # calculate this
    S = Q_(0.5, "dimensionless")
    c = Q_(
        1, "dimensionless"
    )  # the first fraction on pg 2-17 -- essentially the conversion factor
    signal_denom = G_R * C_t * sqrt(power) * B_m * n_B * S * (S + 1) * Q
    signal_denom = signal_denom.to(Q_("G") * sqrt(Q_("W")) * Q_("s"))
    # }}}
    fl.text(strm(
        f"$G_R={G_R:~L}$\n",
        f"$C_t={C_t:~L}$\n",
        f"$power={power:~L}$\n",
        f"$B_m={B_m:~L}$\n",
        f"$Q={Q:~L} $\n",
        f"$n_B={n_B:~L} $\n",
        f"$S={S:~L} $\n",
        f"$c={c:~L} $\n",
        f"$signal\;denom={signal_denom:~L}$",
    ))
    # {{{ configure all the plots -- I like to do this in one place so I
    #     can easily control placement in the PDF
    fl.text(r"\textbf{\texttt{%s}}\par"%file_name)
    fl.next(f"{which_plot} show background subtraction in abs mode")
    gca().set_title("abs mode:\nspectrum and background")
    fl.text(r'\par')
    fl.next(f"{which_plot} absorption, bg. no bl.")
    gca().set_title("abs mode:\nbackground subtracted, show baseline")
    fl.text(r'\par')
    fl.next(f"{which_plot} baseline diagnostic")
    gca().set_title("zoomed-in baseline diagnostic\nshowing only baseline\n(data and fit)")
    # }}}
    d.set_plot_color(next(colors)["color"])
    fl.next(f"{which_plot} show background subtraction in abs mode")
    forplot = background.C.integrate(fieldaxis, cumulative=True)
    forplot.set_plot_color(d.get_plot_color())
    fl.plot(forplot, alpha=0.2)
    forplot = d.C.integrate(fieldaxis, cumulative=True)
    fl.plot(forplot, alpha=0.5, label=label)
    fl.next(f"{which_plot} absorption, bg. no bl.")
    # show the same thing on the main plot as a faint line
    fl.plot(forplot, alpha=0.2, label=label)
    d -= background # subtract background in derivative mode
    d_abs = d.C.integrate(fieldaxis, cumulative=True)
    fl.plot(d_abs, alpha=0.5, label=label)
    # for a protein, we assume we have 1 broad peak that's always over threshold*max
    right_side = abs(d_abs.data)[-1]
    # define a peak as anything that rises above threshold*max
    peaklist = d_abs.contiguous(
        lambda x: abs(x) > right_side + (abs(x).data.max() - right_side) * threshold
    )
    # pull the furthest left and furthest right boundaries of any peaks that we find
    specrange = (peaklist.ravel().min(), peaklist.ravel().max())
    fl.next(f"{which_plot} absorption, bg. no bl.")
    generous_limits = specrange + np.diff(specrange).item() * r_[-pushout, +pushout]
    for j in r_[np.array(specrange)]:
        # show where the peaks were identified
        axvline(x=j, alpha=0.5, color=d.get_plot_color(), ls=":")
    for j in r_[generous_limits]:
        # show the actual limits
        axvline(x=j, alpha=0.5, color=d.get_plot_color(), ls="-")
    d_baseline = d_abs[
        fieldaxis,
        lambda x: np.logical_or(x < generous_limits[0], x > generous_limits[1]),
    ]
    fl.next(f"{which_plot} baseline diagnostic")
    fl.plot(d_baseline, ".", alpha=0.3, label=label, human_units=False)
    # {{{ we don't have a spline in pyspecdata yet, so just hack it
    spl = UnivariateSpline(d_baseline.getaxis(fieldaxis), d_baseline.data)
    polybaseline = d.copy(data=False)
    polybaseline.data = spl(d.getaxis(fieldaxis))
    # }}}
    fl.plot(polybaseline, alpha=0.5, human_units=False)
    fl.next(f"{which_plot} absorption, bg. no bl.")
    fl.plot(polybaseline, alpha=0.5)
    d_abs -= polybaseline
    d_abs.integrate(fieldaxis, cumulative=True)
    fl.next("dblint ÷ denom * conversion")
    d_abs /= signal_denom.magnitude
    d_abs *= dint_conversion
    d_abs /= 1e-6  # convert to μM
    final_conc = (
        d_abs[fieldaxis : (generous_limits[-1], None)].mean(fieldaxis).item()
    ).real
    #print("AVerage concentration:", final_conc)
    d_abs.name("conc").set_units("μM")
    fl.plot(d_abs, alpha=0.5, label=f"{label}, %0.4f μM" % final_conc)
    fl.grid()
    vars[label] = final_conc
    with open(pickle_file, "wb") as fp:
        vars = pickle.dump(vars, fp)
