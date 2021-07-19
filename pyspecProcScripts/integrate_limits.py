from pyspecdata import *
from pylab import subplots, axvline
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like


def integrate_limits(s, axis="t2", fwhm=100, fl=None):
    temp = s.C.mean_all_but(axis)
    if fl is not None:
        fl.next('integration diagnostic')
        fl.push_marker()
        #fl.plot(abs(temp), label='abs before')
    temp_copy = temp.C
    fl.next('before convolution')
    fl.plot(temp_copy)
    if fl is not None:
        fl.plot(temp)

    # pulled from apodization code
    sigma = nddata(np.linspace(1e-10,1e-1,1000),'sigma').set_units('sigma','s')
    gaussians = np.exp(-temp.C.fromaxis(axis)**2/2/sigma**2)
    signal_E = (abs(temp * gaussians)**2).sum(axis)
    signal_E /= signal_E.data.max()
    filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    if fl is not None:
        fl.next('integration diagnostic -- signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    print("*** *** ***")
    print("FILTER WIDTH IS",filter_width)
    print("*** *** ***")
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    # COPY EXACTLY FROM PYSPECDATA CONVOLVE.PY
    for_manual = temp.C
    rough_center = abs(for_manual).C.mean_all_but('t2').argmax('t2').item()
    for_manual.setaxis('t2', lambda t: t- rough_center)
    for_manual.register_axis({'t2':0})
    x = for_manual.C.fromaxis('t2')
    convfunc = lambda x,y: exp(-(x**2)/(2.0*(y**2)))
    myfilter = convfunc(x,filter_width)
    fl.next('Gaussian filter')
    fl.plot(myfilter)
    fl.next('Compare Gaussian Filter and Abs Data Normalized')
    fl.plot(abs(for_manual)/abs(for_manual).max(), alpha=0.5, label='abs data')
    fl.plot(myfilter, alpha=0.5, label='Gaussian filter')
    newdata = for_manual*myfilter
    fl.next('Overlay abs data')
    fl.plot(abs(for_manual), alpha=0.5, label='before applying filter')
    fl.plot(abs(newdata), alpha=0.5, label='after applying filter')
    # END COPY FROM PYSPECDATA CONVOLVE.PY
    #temp.convolve('t2', filter_width)
    if fl is not None:
        fl.next('integration diagnostic')
        fl.plot(abs(temp)/abs(temp).max(), label='after')
        fl.show();quit()
    freq_limits = temp.contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
    if fl is not None:
        fl.next('integration diagnostic')
        axvline(x = freq_limits[0], c='k', alpha=0.75)
        axvline(x = freq_limits[-1], c='k', alpha=0.75)
        fl.pop_marker()
    freq_limits = np.array(freq_limits)
    return freq_limits
