from pylab import*
from pyspecdata import *
from scipy.optimize import leastsq, minimize, basinhopping
import numpy as np
def correl_align(s, align_phases=False,indirect_dim='indirect',fl=None):
    sig_energy = (abs(s)**2).data.mean().item()
    print("*** *** ***")
    print("Before",sig_energy)
    print("*** *** ***")
    if fl is not None:
        fl.next('before correlation\nsig. energy=%g'%sig_energy)
        fl.image(s)
    if align_phases:
        ph0 = s.C.sum('t2')
        ph0 /= abs(ph0)
        s /= ph0
    s.ift('t2')
    correl = s * 0
    for q in range(1,ndshape(s)['vd']):
        temp = s.C
        temp = nddata(roll(temp.data.conj(),q,axis=0),['vd','t2'])
        correl +=s * temp
    if fl is not None:
        fl.next('Look at correlation function - time domain')
        fl.image(correl)
    correl.setaxis('t2',s.getaxis('t2'))
    correl.setaxis('vd',s.getaxis('vd'))
    correl.ft('t2', pad=2**14)
    f_shift = correl.run(real).argmax('t2')
    f_shift = nddata(f_shift.data,[indirect_dim])
    s *= np.exp(-1j*2*pi*f_shift*s.fromaxis('t2'))
    s.ft('t2')
    if fl is not None:
        fl.next('after correlation\nbefore ph0 restore')
        fl.image(s)
    if align_phases:
        s *= ph0
    sig_energy = (abs(s)**2).data.mean().item()
    print("*** *** ***")
    print("after",sig_energy)
    print("*** *** ***")
    if fl is not None:
        fl.next('After correlation\nph0 restored sig energy=%g'%sig_energy)
        fl.image(s)
    return s    
    # {{{ Apodization
    # from Cavanagh p 114
    Kaiser = False
    Gaussian = False
    Hamming = False
    if Kaiser:
        theta = 1.5*pi
        Kaiser_func = special.j0(theta*sqrt(1-(correl.getaxis('t2')**2)/correl.getaxis('t2').max()**2))/special.j0(theta)
        Kaiser_func = nddata(Kaiser_func,'t2')
        Kaiser_func.setaxis('t2',correl.getaxis('t2'))
        if fl is not None:
            fl.next('Look at Kaiser')
            fl.plot(Kaiser_func)
        correl *= Kaiser_func
    if Gaussian:
        mu = 0
        sigma = 0.009
        Gauss_func = exp(-(correl.fromaxis('t2')-mu)**2/(2*sigma**2))
        Gauss_func.setaxis('t2',correl.getaxis('t2'))
        if fl is not None:
            fl.next('Look at Gaussian')
            fl.plot(Gauss_func)
        correl *= Gauss_func
    if Hamming:
        Hamming_func = cos(pi*correl.fromaxis('t2')/(2*correl.getaxis('t2').max()))
        Hamming_func.setaxis('t2',correl.getaxis('t2'))
        if fl is not None:
            fl.next('Look at Hamming')
            fl.plot(Hamming_func)
        correl *= Hamming_func
    #}}}
    if fl is not None:
        fl.next('correlation function\ntime domain, after apod')
        fl.image(correl)
    correl.ft_clear_startpoints('t2')
    correl.ft('t2', shift=True, pad=2**14)
    #correl.ft('t2', pad=2**14)
    if fl is not None:
        fl.next('correlation function\nfreq domain, after apod')
        fl.image(correl)
    f_shift = correl.run(real).argmax('t2')
    f_shift = nddata(f_shift.data,[indirect_dim])
    s *= exp(-1j*2*pi*f_shift*s.fromaxis('t2'))
    s.ft('t2')
    if fl is not None:
        fl.next('after correlation\nbefore ph0 restore')
        fl.image(s)
    if align_phases:
        s *= ph0
    sig_energy = (abs(s)**2).data.mean().item()
    print("*** *** ***")
    print("After",sig_energy)
    print("*** *** ***")
    if fl is not None:
        fl.next('after correlation\nph0 restored sig. energy=%g'%sig_energy)
        fl.image(s)
    return s

