from pyspecdata import *
from pylab import subplots
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging


def integrate_limits(s, axis="t2", fwhm=100, fl=None):
    fwhm = fwhm_calculator(s, axis, fl=fl)
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    d = s.C
    d.mean_all_but(axis).convolve(axis, fwhm/(2*np.sqrt(np.log(4))))
    if fl is not None:
        fl.plot(d/abs(d.data).max(), ax=ax2)
        #fl.pop_marker()
    print(d.getaxis('t2'))
    return d.contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
