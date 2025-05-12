"""
Weighted Integral Function Example
==================================

The weighted integral function is utilized in integrating a signal.
Specifically here, we generate a sinc function whose FT is a heaviside hat
function with a width equal to the integration bounds.
"""
# TODO ☐: pulling very heavily from the stuff in your paper, can we
#         please just show a repeated series of simulated data, and
#         integrate both in the frequency domain and in the time domain.  Do
#         this without apo, and use the same (manual) frequency range for both,
#         so that the results match.

from pyspecdata import r_, nddata, figlist_var, ndshape
from pyspecProcScripts import apod_matched_filter, Heaviside_time_domain
import numpy as np
from pylab import rcParams
from numpy.random import seed
import matplotlib.pyplot as plt

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1

n_repeats = 50
npts = 8192
t = nddata(r_[-npts:npts] * 10 / npts, "t2").set_units("t2", "s")
FWHM = 2 / np.pi
# {{{ Generate data
data = np.pi * FWHM / 2 * np.exp(1j * 2 * np.pi * t - np.pi * FWHM * abs(t))
data = (
    ndshape([len(data["t2"]), n_repeats], ["t2", "repeats"])
    .alloc()
    .setaxis("t2", data.getaxis("t2"))
    .set_units("t2", "s")
    .setaxis("repeats", r_[0:n_repeats])
)
for j in range(len(data.getaxis("repeats"))):
    data["repeats", j].add_noise(0.2)
data.set_plot_color("black")
# }}}
with figlist_var() as fl:
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
    # {{{ Determine integration bounds
    _, filter_width = apod_matched_filter(
        data.C.mean("repeats"), ret_width=True, axis="t2"
    )
    # {{{ make weighted integral function
    filter_width = 4 / filter_width
    f_range = tuple(r_[-filter_width / 2, filter_width / 2])
    wt_fn = Heaviside_time_domain(data, f_range)
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
