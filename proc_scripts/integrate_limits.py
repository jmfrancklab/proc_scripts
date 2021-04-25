from pyspecdata import *
from pylab import subplots
import numpy as np
import logging


def integrate_limits(s, axis="t2", fwhm=100, fl=None):
    signal_sign = s.C.sum(axis).run(np.real).run(np.sign)
    print(ndshape(s))
    temp = s.real * signal_sign
    #fl.push_marker()
    if fl is not None:
        fig, (ax1,ax2) = subplots(2,1)
        fl.next("integration diagnostic", fig=fig)
        fl.plot(temp, ax=ax1)
    temp.mean_all_but(axis)
    if fl is not None:
        fl.plot(temp/abs(temp.data).max(), ax=ax2)
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    temp.convolve(axis, fwhm/(2*np.sqrt(np.log(4))))
    if fl is not None:
        fl.plot(temp/abs(temp.data).max(), ax=ax2)
    #fl.pop_marker()
    return temp.contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
