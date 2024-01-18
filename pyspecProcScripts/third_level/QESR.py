import numpy as np
from numpy import r_
import os
from matplotlib.pyplot import axvline, gca
import matplotlib.pyplot as plt
from pyspecdata import find_file, gammabar_e, ndshape
from scipy.interpolate import UnivariateSpline
from ..first_level.QESR_rescale import QESR_scalefactor
import pickle

# {{{ all basic info
colors = plt.rcParams[
    "axes.prop_cycle"
]()  # this is the default matplotlib cycler for line styles
fieldaxis = "$B_0$"
# }}}

"""Process an ESR spectra
=========================
Converts a raw ESR spectra to the double integral and calculates the
appropriate concentration of spins present. 

If a background spectra is provided, background subtraction is applied to
the loaded ESR spectra. The double integral is finally rescaled and
multiplied by a proportionality constant associated with the matching value
in your pyspecdata config file.  The resulting plot shows the double
integral with the calculated concentration in micromolar in the legend.
"""

def QESR(file_name, 
        label, 
        pushout=0.5,
        threshold=0.05, 
        pickle_file=None, 
        background=None,
        exp_type=None,
        which_plot=None, 
        calibration_name=None,
        diameter_name=None, 
        color=None,
        fl=None):
    """
    Parameters
    ==========
    pushout: float
        Adjusts the "generous limits" factor which is added
        to the integration area surrounding the peaks
    threshold: float
        Used in defining the integration limits (threshold*max)
    pickle_file: str
        Name of pickle file you want to append the concentration to
    background: nddata
        The background spectrum that should already be loaded
    exp_type: str
        Where the file of interest is
    which_plot: str
        Name that will be used in titles of plots. This is useful when you
        are plotting more than 1 dataset
    calibration_name:  str
        Name of the value in your pyspecdata config file that points to
        the proportionality constant
    diameter_name:  str
        Name of the value in your pyspecdata config file that points to 
        the diameter of the capillary being used
    color:  str
            String of the desired color for the plot. This is especially
            useful when trying to coordinate color matching with other
            plots
    fl: figure list
        required!
    """
    if exp_type is None:
        raise ValueError("You must specify the location of the file with exp_type!")
    if which_plot is None:
        which_plot = file_name
    if fl is None:
        raise ValueError("for now, you just have to pass a figure list")
    if pickle_file is not None:
        # {{{ make a pickle file, that contains a
        #     dictionary of the concentrations of all samples
        #     that we are looking at
        #     the label should have the batch date and
        #     residue, so should be good identifier
        if os.path.exists(pickle_file):
            with open(pickle_file, "rb") as fp:
                pickle_vars = pickle.load(fp)
        else:
            pickle_vars = {}
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
    # {{{ configure all the plots -- I like to do this in one place so I
    #     can easily control placement in the PDF
    fl.text(r"\textbf{\texttt{%s}}\par"%file_name)
    fl.next(f"{which_plot} show background subtraction in abs mode")
    gca().set_title("abs mode:\nspectrum and background spectrum")
    fl.text(r'\par')
    fl.next(f"{which_plot} absorption, bg. no bl.")
    gca().set_title("abs mode:\nbackground subtracted, show baseline")
    fl.text(r'\par')
    fl.next(f"{which_plot} baseline diagnostic")
    gca().set_title("zoomed-in baseline diagnostic\nshowing only baseline\n(data and fit)")
    # }}}
    if color is None:
        color = next(colors)["color"]
    d.set_plot_color(color)
    fl.next(f"{which_plot} show background subtraction in abs mode")
    forplot = background.C.integrate(fieldaxis, cumulative=True)
    forplot.set_plot_color(d.get_plot_color())
    fl.plot(forplot, alpha=0.2)
    forplot = d.C.integrate(fieldaxis, cumulative=True)
    fl.plot(forplot, alpha=0.5, label=label)
    fl.next(f"{which_plot} absorption, bg. no bl.")
    # show the same thing on the main plot as a faint line
    fl.plot(forplot, alpha=0.2, label="%s - no background spec. subtracted"%label)
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
    spl = UnivariateSpline(d_baseline.getaxis(fieldaxis), d_baseline.data.real)
    polybaseline = d.copy(data=False)
    polybaseline.data = spl(d.getaxis(fieldaxis))
    # }}}
    fl.plot(polybaseline, alpha=0.5, human_units=False)
    fl.next(f"{which_plot} absorption, bg. no bl.")
    fl.plot(polybaseline, alpha=0.5, label="spline/poly baseline\n(integrate between this and other curve)", lw=1)
    d_abs -= polybaseline #background subtraction
    d_abs.integrate(fieldaxis, cumulative=True)
    fl.next("dblint ÷ denom * conversion")
    d_abs /= QESR_scalefactor(d, calibration_name=calibration_name,
            diameter_name=diameter_name)
    final_conc = (
        d_abs[fieldaxis : (generous_limits[-1], None)].mean(fieldaxis).item()
    ).real
    d_abs.name("conc").set_units("μM")
    fl.plot(d_abs, alpha=0.5, label=f"{label}, %0.4f μM" % final_conc)
    fl.grid()
    if pickle_file is not None:
        pickle_vars[label] = final_conc
        with open(pickle_file, "wb") as fp:
            pickle_vars = pickle.dump(pickle_vars, fp)
    return final_conc
