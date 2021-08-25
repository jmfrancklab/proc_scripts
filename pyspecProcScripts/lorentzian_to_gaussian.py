from pyspecdata import *
from .apod_matched_filter import apod_matched_filter

def L2G(s, axis='t2',
        fl=None):
    r"""
    Lorentzian-to-Gaussian
    ============================

    This performs Lorentzian-to-Gaussian transformation on input data. Data
    should be input as time domain data.

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
    _,filter_width = apod_matched_filter(temp,
            convolve_method='lorentzian',
            ret_width=True,
            basename='L2G',
            fl=fl)
    if fl is not None:
        fl.next('Lorentzian to Gaussian diagnostic -- time domain')
        fl.plot(temp,
                alpha=0.6,
                label='before')
        fl.next('Lorentzian to Gaussian diagnostic -- frequency domain')
        fl.plot(temp.C.ft('t2'),
                alpha=0.6,
                label='before')
    temp *= np.exp(-temp.fromaxis(axis)**2/2/(filter_width*pi/2)**2)
    temp /= np.exp(-abs(temp.fromaxis(axis))/filter_width)
    if fl is not None:
        fl.next('Lorentzian to Gaussian diagnostic -- time domain')
        fl.plot(temp,
                alpha=0.6,
                label='after')
        fl.next('Lorentzian to Gaussian diagnostic -- frequency domain')
        fl.plot(temp.C.ft('t2'),
                alpha=0.6,
                label='after')
    return temp
