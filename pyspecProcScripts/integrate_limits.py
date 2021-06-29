from pyspecdata import *
from pylab import subplots, axvline
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging


def integrate_limits(s, axis="t2", fwhm=100, fl=None):
    temp = s.C.mean_all_but(axis)
    if fl is not None:
        fl.next('integration diagnostic')
        fl.push_marker()
        fl.plot(abs(temp))
    temp.run(np.real)
    temp /= temp.C.sum(axis).run(np.sign) # flip it so it's pointing up
    if fl is not None:
        fl.plot(temp)

    # pulled from apodization code
    sigma = nddata(np.linspace(1e-5,1e3,1000),'sigma').set_units('sigma','s')
    gaussians = np.exp(-temp.C.fromaxis(axis)**2/2/sigma**2)
    signal_E = (abs(temp * gaussians)**2).sum(axis)
    signal_E /= signal_E.data.max()
    filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    if fl is not None:
        fl.next('integration diagnostic -- signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    fwhm = filter_width
    print("FWHM IS",fwhm)
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    temp.convolve(axis, fwhm/(2*np.sqrt(np.log(4))))
    #temp *= sqrt(fwhm/(2*np.sqrt(np.log(4))))
    if fl is not None:
        fl.next('integration diagnostic')
        fl.plot(temp)
    freq_limits = temp.contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
    if fl is not None:
        fl.next('integration diagnostic')
        axvline(x = freq_limits[0], c='k', alpha=0.75)
        axvline(x = freq_limits[-1], c='k', alpha=0.75)
        fl.pop_marker()
    return freq_limits
