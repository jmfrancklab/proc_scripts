from pyspecdata import *
def align_and_slice(s, axis='t2', convwidth=500, threshold=0.05, fl=None):
    r"""Align the peak frequencies, and automatically slice them -- passing fl assumes you are debugging and want diagnostic plots"""
    if fl is not None:
        fl.push_marker()
    if not s.get_ft_prop(axis):
        raise ValueError("Data has not been FTd yet!")
    s.ift(axis)
    R = pi*convwidth # Lorentzian FWHM formula -- R = 1/T2, FWHM = R/pi
    s_conv = s.C*exp(-R * abs(s.fromaxis(axis)))
    s_conv_save = s.C
    s_conv.ft(axis)
    if fl is not None:
        fl.next('check for center')
        fl.image(abs(s_conv))
    center_frq = abs(s_conv).argmax(axis)
    s *= exp(-1j*2*pi*center_frq*s.fromaxis(axis))
    s_conv = s_conv_save * exp(-1j*2*pi*center_frq*s.fromaxis(axis))
    s.ft(axis)
    s_conv.ft(axis)
    s_conv.mean_all_but([axis])
    slices = s_conv.contiguous(lambda x:
            abs(x)>threshold*abs(s_conv.data).max())
    if fl is None:
        assert slices.shape[0] == 1, "found more than one peak: "+str(slices)
    slices = tuple(slices[0,:])
    print("slicing to",slices)
    if fl is not None:
        fl.pop_marker()
    return s[axis:slices]
def correlation_align(s,avg,convwidth=0,axis='t2',fl=None):
    assert not s.get_ft_prop(axis), 'I want time-domain data'
    assert not avg.get_ft_prop(axis), 'I want time-domain data'
    avg.run(conj)
    R = pi*convwidth # Lorentzian FWHM formula -- R = 1/T2, FWHM = R/pi
    # ALSO convolve both
    check = avg * s * exp(-2*R*abs(s.fromaxis(axis)))
    assert s.get_ft_prop(axis,['start','freq']) is not None, "should be FT'd first, so the startpoint is set"
    check.ft(axis, pad=2**14)
    check.reorder(axis,first=False)
    forplot = abs(check)[axis:(-50,50)]
    myline = forplot.C.argmax(axis)
    indirect_dims = set(s.dimlabels) - set([axis])
    phcyc_dims = [j for j in indirect_dims if j.startswith('ph')]
    phcyc_dims.sort()
    indirect_dims = list(set(indirect_dims) - set(phcyc_dims))
    #forplot.reorder(phcyc_dims+indirect_dims)
    forplot.smoosh(phcyc_dims+indirect_dims,'indirect'
            ).setaxis('indirect','#'
                    ).reorder('indirect',first=False)
    if fl is not None:
        fl.push_marker()
        fl.next('cross-correlation')
        fl.image(forplot,human_units=False)
        fl.plot(myline.smoosh(phcyc_dims+indirect_dims,'indirect'
            ).setaxis('indirect','#'
                    ).reorder('indirect',first=True)
                , 'w', linewidth=3, alpha=0.5,
                human_units=False)
    return s*exp(-1j*2*pi*myline*s.fromaxis(axis))
