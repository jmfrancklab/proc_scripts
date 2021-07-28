from pyspecdata import *
from pylab import subplots, axvline
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like
from matplotlib.pyplot import annotate


def integrate_limits(s, axis="t2", fwhm=100, fl=None):
    temp = s.C.mean_all_but(axis)
    if fl is not None:
        fl.next('integration diagnostic')
        fl.push_marker()
        fl.plot(abs(temp.C.ft('t2'))/abs(temp.C.ft('t2')).max(), alpha=0.6, label='before Lorentzian to Gauss')
    temp_copy = temp.C
    fl.next('before convolution')
    fl.plot(temp_copy)
    if fl is not None:
        fl.plot(temp)

    # pulled from apodization code
    Gaussian_Conv = False
    Lorentzian_Conv = True
    sigma = nddata(np.linspace(1e-10,1e-1,1000),'sigma').set_units('sigma','s')
    if Gaussian_Conv:
        gaussians = np.exp(-temp.C.fromaxis(axis)**2/2/sigma**2)
        signal_E = (abs(temp * gaussians)**2).sum(axis)
    if Lorentzian_Conv:
        lorentzians = np.exp(-abs(temp.C.fromaxis(axis))/sigma)
        signal_E = (abs(temp * lorentzians)**2).sum(axis)
    signal_E /= signal_E.data.max()
    if Gaussian_Conv:
        filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    if Lorentzian_Conv:
        filter_width = abs(signal_E-signal_E.max()/2).argmin('sigma').item()
    if fl is not None:
        fl.next('integration diagnostic -- signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    print("*** *** ***")
    print("FILTER WIDTH IS",filter_width)
    print("*** *** ***")
    if Gaussian_Conv:
        convfunc = lambda x,y: exp(-(x**2)/(2.0*(y**2)))
    if Lorentzian_Conv:
        convfunc = lambda x,y: exp(-(abs(x))/y)
    Gaussian_func = lambda x,y: exp(-(x**2)/(2.0*(y**2)))
    Lorentzian_func = lambda x,y: exp(-(abs(x))/y)
    # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    # COPY EXACTLY FROM PYSPECDATA CONVOLVE.PY
    rough_center = abs(temp).C.mean_all_but('t2').argmax('t2').item()
    temp.setaxis('t2', lambda t: t- rough_center).register_axis({'t2':0})
    #temp = temp['t2':(0.5e-3,None)]
    temp = temp['t2':(7e-3,None)]
    manual_convolve = False
    if manual_convolve:
        # assumed temp starts in time domain
        for_manual = temp.C
        x = for_manual.C.fromaxis('t2')
        myfilter = convfunc(x,filter_width)
        if Gaussian_Conv:
            fl.next('Gaussian filter')
        if Lorentzian_Conv:
            fl.next('Gaussian filter')
        fl.plot(myfilter)
        fl.show();quit()
        if Gaussian_Conv:
            fl.next('Compare Gaussian Filter and Abs Data Normalized')
        if Lorentzian_Conv:
            fl.next('Compare Lorentzian Filter and Abs Data Normalized')
        fl.plot(abs(for_manual)/abs(for_manual).max(), alpha=0.5, label='abs data')
        fl.plot(myfilter, alpha=0.5, label='filter')
        newdata = for_manual*myfilter
        fl.next('Overlay abs data')
        fl.plot(abs(for_manual), alpha=0.5, label='before applying filter')
        fl.plot(abs(newdata), alpha=0.5, label='after applying filter')
        temp = newdata.C
        temp.ft('t2')
    # END COPY FROM PYSPECDATA CONVOLVE.PY
    if not manual_convolve:
        temp.ft('t2')
        print("I want this filter_width",filter_width)
        # filter_width is λ/π which is FWHM for Lorentzian in Hz
        # filter for L-to-G lies between λ/2 and λ*2 (Cav p.116)
        # thus we need to multiply our filter by π, in order to have filters of
        # the appropriate range
        this_filter = (filter_width*pi)/2
        print("Filter width for Lorentz-to-Gauss",this_filter)
        temp.convolve('t2', this_filter, convfunc=Gaussian_func)
        temp.ift('t2')
        temp /= Lorentzian_func(temp.C.fromaxis('t2'),filter_width)
        temp.ft('t2')
        #temp.convolve('t2', filter_width, convfunc=convfunc)
    if fl is not None:
        fl.next('integration diagnostic')
        fl.plot(abs(temp)/abs(temp).max(), alpha=0.6, label='after Lorentzian-to-Gauss')
        fl.show();quit()
    limit_for_contiguous = 0.125
    freq_limits = temp.contiguous(lambda x: abs(x) > limit_for_contiguous * abs(x).data.max())[0]
    if fl is not None:
        fl.next('integration diagnostic')
        axvline(x = freq_limits[0], c='k', alpha=0.75)
        axvline(x = freq_limits[-1], c='k', alpha=0.75)
        annotate(str(limit_for_contiguous), xy=(freq_limits[-1],0.85))
        fl.pop_marker()
    freq_limits = np.array(freq_limits)
    return freq_limits
