from pyspecdata import *
from pylab import subplots
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging


def integrate_limits(s, axis="t2", fwhm=100, fl=None):
    signal_sign = s.C.sum(axis).run(np.real).run(np.sign)
    temp = s.real * signal_sign

    # pulled from apodization code
    sigma = nddata(np.linspace(1e-5,1e3,1000),'sigma').set_units('sigma','s')
    s_avg = s.C.mean_all_but('t2')
    gaussians = np.exp(-s_avg.C.fromaxis('t2')**2/2/sigma**2)
    signal_E = (abs(s_avg * gaussians)**2).sum('t2')
    signal_E /= signal_E.data.max()
    filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    if fl is not None:
        fl.next('signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    fwhm = filter_width
    fwhm += 50
    print("FWHM IS",fwhm)
    fl.push_marker()
    if fl is not None:
        fig, (ax1,ax2) = subplots(2,1)
        fl.next("integration diagnostic", fig=fig)
        fl.plot(temp, ax=ax1)
    temp.mean_all_but(axis)
    if fl is not None:
        fl.plot(temp/abs(temp.data).max(), ax=ax2,human_units=False)
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    temp.convolve(axis, fwhm/(2*np.sqrt(np.log(4))))
    if fl is not None:
        fl.plot(temp/abs(temp.data).max(), ax=ax2,human_units=False)
    fl.pop_marker()
    return temp.contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
