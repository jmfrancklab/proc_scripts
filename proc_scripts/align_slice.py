from pyspecdata import *
import matplotlib.pyplot as plt
import numpy as np
def align_and_slice(s, dimension='t2', convwidth=500, threshold=0.05, fl=None):
    r"""Align the peak frequencies, and automatically slice them -- passing fl assumes you are debugging and want diagnostic plots"""
    if fl is not None:
        fl.push_marker()
    if not s.get_ft_prop(dimension):
        raise ValueError("Data has not been FTd yet!")
    s.ift(dimension)
    R = pi*convwidth # Lorentzian FWHM formula -- R = 1/T2, FWHM = R/pi
    s_conv = s.C*exp(-R * abs(s.fromaxis(dimension)))
    s_conv_save = s.C
    s_conv.ft(dimension)
    if fl is not None:
        fl.next('check for center')
        fl.image(abs(s_conv))
    center_frq = abs(s_conv).argmax(dimension)
    s *= exp(-1j*2*pi*center_frq*s.fromaxis(dimension))
    s_conv = s_conv_save * exp(-1j*2*pi*center_frq*s.fromaxis(dimension))
    s.ft(dimension)
    s_conv.ft(dimension)
    s_conv.mean_all_but([dimension])
    slices = s_conv.contiguous(lambda x:
            abs(x)>threshold*abs(s_conv.data).max())
    if fl is None:
        assert slices.shape[0] == 1, "found more than one peak: "+str(slices)
    slices = tuple(slices[0,:])
    print("slicing to",slices)
    if fl is not None:
        fl.pop_marker()
    return s[axis:slices]
def correlation_align(s,avg,convwidth=0,axis='t2',color='k',linestyle='',fl=None):
    assert not s.get_ft_prop(axis), 'I want time-domain data'
    assert not avg.get_ft_prop(axis), 'I want time-domain data'
    avg.run(conj)
    R = pi*convwidth # Lorentzian FWHM formula -- R = 1/T2, FWHM = R/pi
    # ALSO convolve both
    check = avg * s * exp(-2*R*abs(s.fromaxis(axis)))
    assert s.get_ft_prop(axis,['start','freq']) is not None, "should be FT'd first, so the startpoint is set"
    check.ft(axis, pad=2**14)
    #check.reorder(axis,first=False) #does nothing, as in with this statement uncommented looks the same as 
    #with it commented out
    forplot = abs(check)[axis:(-500,500)]
    myline = forplot.C.argmax(axis)
    indirect_dims = set(s.dimlabels) - set([axis])
    phcyc_dims = [j for j in indirect_dims if j.startswith('ph')]
    phcyc_dims.sort()
    indirect_dims = list(set(indirect_dims) - set(phcyc_dims))
    forplot.smoosh(phcyc_dims+indirect_dims,'indirect',noaxis=True).setaxis('indirect','#').reorder('indirect',first=False)
    thisline = myline.C
    thisline.smoosh(phcyc_dims+indirect_dims, 'indirect', noaxis=True).setaxis('indirect','#').reorder('indirect',first=True)
    if fl is not None:
        fl.push_marker()
        fl.next('cross-correlation')
        fl.image(forplot,human_units=False)
        fl.plot(thisline.imag,color=color,linestyle='-', linewidth=3, alpha=0.5,
                human_units=False)
        fl.plot(thisline.real, color='k', linestyle='--', linewidth=3, alpha=0.5,
                human_units=False)
    if fl is not None:
        fl.next('s and avg')
        fl.plot(s.real,color='r',human_units=False)
        fl.plot(s.imag,color='b',human_units=False)
        fl.plot(avg.real,color='m',human_units=False)
        fl.plot(avg.imag,color='g', human_units=False)
    if fl is not None:
        fl.next('cross-correlation plot')
        fl.plot(thisline.real, color='k', linestyle='--', linewidth=3, alpha=0.5,
                human_units=False)
        fl.plot(thisline.imag, color='r', linestyle='-', linewidth=3, alpha=0.5,
                human_units=False)
    return s*exp(-1j*2*pi*myline*s.fromaxis(axis))

