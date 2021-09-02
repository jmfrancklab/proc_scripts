from pyspecdata import *
from pylab import subplots, axvline
import numpy as np
from .fwhm_calculate import fwhm_calculator
import logging
from pylab import r_,fft,ifft,ifftshift,fftshift,exp,ones_like
from matplotlib.pyplot import annotate
from .apod_matched_filter import apod_matched_filter

def integrate_limits(s, axis="t2",
        convolve_method='gaussian',
        cutoff = 0.25,
        fl=None):
    r"""This function takes data in the frequency
    domain and finds the corresponding frequency
    limits of the signal for integration.

    Limits are found by applying a matched filter
    (choice of Lorentzian or Gaussian ``convolve_method``) to the
    frequency domain data, which it then uses to
    determine the integration limits based on a
    ``cutoff`` amplitude based as a fraction of the maximum signal intensity.

    Parameters
    ==========

    axis: str
        apply convolution along `axisname`

    convolve_method: str
        specify as one of the following 3 options:
        * Option 1 - 'gaussian' 
        to calculate the matched-filter via
        Gaussian convolution of the frequency
        domain signal.
        * Option 2 - 'lorentzian' 
        to calculate the matched-filter via
        Lorentzian convolution of the frequency
        domain signal.

    cutoff: float
        set the integration limit to the points where the convolved signal
        crosses this fraction of its maximum.
        (e.g. cutoff=0.5 would give integration limits that correspond to the
        frequencies marking the full width half maximum of the convolution of
        the signal and its matched filter)

    fl: None or figlist_var()
        to show diagnostic plots, set `fl` to the figure list; set `fl` to None in
        order not to see any diagnostic plots

    """
    assert s.get_ft_prop(axis), "data must be in the frequency domain along %s!!"%axis
    temp = s.C.mean_all_but(axis).real
    # {{{ taking real in the frequency domain enforces symmetry in the
    # time domain.  Without the following, the negative time components
    # get aliased to large time.
    # Might make sense to rather have `.real` above do this, since the
    # reasoning should always apply!  (discuss)
    # AB: Yes, then I agree `.real` should be doing this.
    temp.set_ft_prop(axis,['start','time'],-s.get_ft_prop(axis,'dt')*(ndshape(s)[axis]//2))
    # }}}
    convolve_method = convolve_method.lower()
    if fl is not None:
        fl.push_marker()
        forplot = temp
        fl.next('integration diagnostic')
        fl.plot(forplot/forplot.data.max(), alpha=0.6, label='before convolve')
    temp.ift('t2')
    temp = apod_matched_filter(temp,
            convolve_method=convolve_method,
            fl=fl)
    temp.ft('t2')
    if fl is not None:
        fl.next('integration diagnostic')
        fl.plot(temp/temp.max(), alpha=0.6, label='after convolve')
    freq_limits = temp.contiguous(lambda x: x.real > cutoff * x.real.data.max())[0]
    if fl is not None:
        fl.next('integration diagnostic')
        axvline(x = freq_limits[0], c='k', alpha=0.75)
        axvline(x = freq_limits[-1], c='k', alpha=0.75)
        annotate(str(limit_for_contiguous), xy=(freq_limits[-1],0.85))
        fl.pop_marker()
    freq_limits = np.array(freq_limits)
    # Think still need to return 'temp' if lorentzian_to_gaussian 
    return freq_limits
