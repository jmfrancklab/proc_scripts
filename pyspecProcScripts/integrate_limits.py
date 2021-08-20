from pyspecdata import *
from pylab import subplots, axvline
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like
from matplotlib.pyplot import annotate
logger = init_logging("debug")

def integrate_limits(s, axis="t2",
        filter_width=100,
        convolve_method='Gaussian',
        fl=None):
    r"""
    Integrate Limits
    ============================

    This function takes data in the time domain and applies a matched filter
    (choice of Lorentzian or Gaussian, default is Lorentzian) to the frequency
    domain data, which it then uses to determine the integration limits based on a
    cut off from the maximum signal intensity.

    axisname: str
        apply convolution along `axisname`

    fwhm: int
        width of the matched filter is `fwhm` - this is calculated more precisely
        within the program

    convolve_method: str
        specify as one of the following 3 options:
        * Option 1 - 'Gaussian' 
        to calculate the matched-filter via
        Gaussian convolution of the frequency
        domain signal.
        * Option 2 - 'Lorentzian' 
        to calculate the matched-filter via
        Lorentzian convolution of the frequency
        domain signal.
        * Option 3 - 'Lorentzian_to_Gaussian'
        apply Lorentzian to Gaussian
        transformation to your data, in order to
        determine integral limits.

    fl: None or figlist_var()
        to show diagnostic plots, set `fl` to the figure list; set `fl` to None in
        order not to see any diagnostic plots

    """
    Lorentz_to_Gauss = False
    temp = s.C.mean_all_but(axis)
    if fl is not None:
        fl.next('integration diagnostic')
        fl.push_marker()
        fl.plot(abs(temp.C.ft('t2'))/abs(temp.C.ft('t2')).max(), alpha=0.6, label='before convolve')
    sigma = nddata(np.linspace(1e-10,1e-1,1000),'sigma').set_units('sigma','s')
    if convolve_method == 'Gaussian':
        convolution_set = np.exp(-temp.C.fromaxis(axis)**2/2/sigma**2)
    elif convolve_method == 'Lorentzian':
        convolution_set = np.exp(-abs(temp.C.fromaxis(axis))/sigma)
    signal_E = (abs(temp * convolution_set)**2).sum(axis)
    signal_E /= signal_E.data.max()
    if convolve_method == 'Gaussian':
        filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    elif convolve_method == 'Lorentzian':
        filter_width = abs(signal_E-signal_E.max()/2).argmin('sigma').item()
    logger.info(strm("FILTER WIDTH IS",filter_width))
    if fl is not None:
        fl.next('integration diagnostic -- signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    if convolve_method == 'Gaussian':
        temp *= np.exp(-temp.C.fromaxis(axis)**2/2/filter_width**2)
    temp.ft('t2')

    if Lorentz_to_Gauss:
        # not immediately sure why, but I need to
        # do this in order to get freq limits
        # around the peak (otherwise they are
        # around 0)
        rough_center = abs(temp).C.mean_all_but('t2').argmax('t2').item()
        temp.setaxis('t2', lambda t: t- rough_center).register_axis({'t2':0})
        temp = temp['t2':(7e-3,None)]
    if Lorentz_to_Gauss:
        # filter_width is λ/π which is FWHM for Lorentzian in Hz
        # filter for L-to-G lies between λ/2 and λ*2 (Cav p.116)
        # thus we need to multiply our filter by π, in order to have filters of
        # the appropriate range
        this_filter = (filter_width*pi)/2
        logger.info(strm("Filter width for Lorentz-to-Gauss",this_filter))
        temp.convolve('t2', 1/this_filter, convfunc=Gaussian_func)
        temp.ift('t2')
        temp /= Lorentzian_func(temp.C.fromaxis('t2'),filter_width)
        temp.ft('t2')
    if fl is not None:
        fl.next('integration diagnostic')
        fl.plot(abs(temp)/abs(temp).max(), alpha=0.6, label='after convolve')
    limit_for_contiguous = 0.25
    freq_limits = temp.contiguous(lambda x: abs(x) > limit_for_contiguous * abs(x).data.max())[0]
    if fl is not None:
        fl.next('integration diagnostic')
        axvline(x = freq_limits[0], c='k', alpha=0.75)
        axvline(x = freq_limits[-1], c='k', alpha=0.75)
        annotate(str(limit_for_contiguous), xy=(freq_limits[-1],0.85))
        fl.pop_marker()
    freq_limits = np.array(freq_limits)
    if Lorentz_to_Gauss:
        return freq_limits,temp
    else:
        return freq_limits
