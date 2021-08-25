from pyspecdata import *

def apod_matched_filter(s, axis='t2',
        convolve_method='gaussian',
        fl=None):
    temp = s.C
    sigma = nddata(np.linspace(1e-10,1e-1,1000),'sigma').set_units('sigma','s')
    if convolve_method == 'gaussian':
        convolution_set = np.exp(-temp.fromaxis(axis)**2/2/sigma**2)
    elif convolve_method == 'lorentzian':
        convolution_set = np.exp(-abs(temp.fromaxis(axis))/sigma)
    signal_E = (abs(temp * convolution_set)**2).sum(axis)
    signal_E /= signal_E.data.max()
    if convolve_method == 'gaussian':
        filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    elif convolve_method in ['lorentzian','lorentzian_to_gaussian']:
        filter_width = abs(signal_E-signal_E.max()/2).argmin('sigma').item()
    logger.info(strm("FILTER WIDTH IS",filter_width))
    if fl is not None:
        fl.next('integration diagnostic -- signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    if fl is not None:
        fl.next('integration diagnostic -- time domain')
        fl.plot(abs(temp), alpha=0.6, label='before mult')
    if convolve_method == 'gaussian':
        temp *= np.exp(-temp.fromaxis(axis)**2/2/filter_width**2)
    elif convolve_method in 'lorentzian':
        temp *= np.exp(-abs(temp.fromaxis(axis))/filter_width)
    return temp
