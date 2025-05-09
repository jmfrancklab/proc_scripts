"""
Weighted Integral Function Example
==================================

The weighted integral function is utilized in integrating a signal.
Specifically here, we generate a sinc function whose FT is a heaviside hat
function with a width equal to the integration bounds.
"""
# TODO ☐: pulling very heavily from the stuff in your paper, can we
#         please just show a repeated series of simulated data, and
#         integrate both in the frequency domain and in the time domain.
#         Do this without apo, and use the same (manual) frequency range for both,
#         so that the results match.

import pyspecdata as psd
from pyspecdata import r_
import pyspecProcScripts as psdpr
import numpy as np
from pylab import rcParams
from numpy.random import seed
import matplotlib.pyplot as plt

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1

npts = 8192
t = psd.nddata(r_[-npts:npts] * 10 / npts, "t2").set_units("t2", "s")
FWHM = 10
f_range = (-50, 50)
# {{{ Generate data
data = np.pi * FWHM / 2 * np.exp(1j * 2 * np.pi * t - np.pi * FWHM * abs(t))
data.set_plot_color("black")
# }}}
with psd.figlist_var() as fl:
    # {{{ Set up plots
    fig, thisax = plt.subplots(2, 1)
    # TODO ☐: this is not a 'weighted integral function'
    fl.next("Weighted Integral Functions", fig=fig)
    fl.skip_units_check()
    thisax[0].set_title("Frequency Domain")
    thisax[1].set_title("Time Domain")
    fig.suptitle("Weighted Integral Functions")
    # }}}
    data.set_units("t2", "s")
    data.ft("t2", shift=True)
    # {{{ make weighted integral function
    wt_fn = psdpr.Heaviside_time_domain(data, f_range)
    wt_fn.set_units("t2", "s")
    wt_fn.ft("t2")
    wt_fn.set_plot_color("tab:blue")
    # }}}
    # {{{ Plot
    for thisdata in [data, wt_fn]:
        fl.plot(thisdata, ax=thisax[0])
        thisdata.ift("t2")
        # Just normalize by the max so they are scaled correctly
        fl.plot(thisdata["t2":(-0.6, 0.6)] / thisdata.max(), ax=thisax[1])
    # }}}
