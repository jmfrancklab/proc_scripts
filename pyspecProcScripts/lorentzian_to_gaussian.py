from pyspecdata import *
from pylab import axhline
from .apod_matched_filter import apod_matched_filter

def L2G(s, axis='t2',
        fl=None):
    r"""
    Lorentzian-to-Gaussian
    ============================

    This performs Lorentzian-to-Gaussian transformation on input data. Data
    should be input as frequency domain data.

    s: nddata
        time domain data to undergo transformation

    axis: str
        transformation along `axisname` (t2)

    fl: None or figlist_var()
        to show diagnostic plots, set `fl` to the
        figure list; set `fl` to None in order not
        to see any diagnostic plots
    """
    temp = s.C.mean_all_but(axis)
    assert temp.get_prop('FT')[axis], "Data must be in frequency domain"
    if fl is not None:
        fl.next('Lorentzian to Gaussian diagnostic -- frequency domain')
        fl.plot(temp,
                alpha=0.6,
                label='before')
    temp.ift(axis)
    _,filter_width = apod_matched_filter(temp,
            convolve_method='lorentzian',
            ret_width=True,
            fl=fl)
    if fl is not None:
        fl.next('Lorentzian to Gaussian diagnostic -- time domain')
        fl.plot(temp,
                alpha=0.6,
                label='before')
    temp *= np.exp(-temp.fromaxis(axis)**2/2/(filter_width*pi/2)**2)
    temp /= np.exp(-abs(temp.fromaxis(axis))/filter_width)
    if fl is not None:
        fl.next('Lorentzian to Gaussian diagnostic -- time domain')
        fl.plot(temp,
                alpha=0.6,
                label='after')
    temp.ft(axis)
    if fl is not None:
        fl.next('Lorentzian to Gaussian diagnostic -- frequency domain')
        fl.plot(temp,
                alpha=0.6,
                label='after')
        axhline(y=0,color='k')
    # {{{ actually apply the transform to the data and return
    temp = s.C.ift(axis)
    temp *= np.exp(-temp.fromaxis(axis)**2/2/(filter_width*pi/2)**2)
    temp /= np.exp(-abs(temp.fromaxis(axis))/filter_width)
    temp.ft(axis)
    # }}}
    return temp
