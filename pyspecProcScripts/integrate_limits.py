from pyspecdata import *
from pylab import subplots
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging


def integrate_limits(s, axis="t2", fl=None):
    """automatically determine the integration limits

    Parameters
    ==========
    axis: str
        name of the dimension along which you want to determine the integration
        limits

    Returns
    =======
    retval: ndarray
        A len 2 ndarray tuple with the start and stop indeces for the integration
        bounds.
    """
    signal_sign = s.C.sum(axis).run(np.real).run(np.sign)
    temp = abs(s).real * signal_sign

    # pulled from apodization code
    sigma = nddata(np.linspace(1e-5,1e3,1000),'sigma').set_units('sigma','s')
    s_avg = s.C.mean_all_but('t2')
    gaussians = np.exp(-s_avg.C.fromaxis('t2')**2/2/sigma**2)
    signal_E = (abs(s_avg * gaussians)**2).sum('t2')
    signal_E /= signal_E.data.max()
    filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    if fl is not None:
        fl.push_marker()
        fl.next('signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
        fl.pop_marker()
    fwhm = filter_width
    fl.push_marker()
    temp.mean_all_but(axis)
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    temp.convolve(axis, fwhm/(2*np.sqrt(np.log(4))))
    fl.pop_marker()
    return temp.contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
