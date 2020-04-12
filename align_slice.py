from pyspecdata import *
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
    return s[dimension:slices]
def correlation_align(s,avg,convwidth=0,dimension='t2',fl=None):
    assert not s.get_ft_prop(dimension), 'I want time-domain data'
    assert not avg.get_ft_prop(dimension), 'I want time-domain data'
    avg.run(conj)
    R = pi*convwidth # Lorentzian FWHM formula -- R = 1/T2, FWHM = R/pi
    # ALSO convolve both
    check = avg * s * exp(-2*R*abs(s.fromaxis(dimension)))
    assert s.get_ft_prop(dimension,['start','freq']) is not None, "should be FT'd first, so the startpoint is set"
    check.ft('t2', pad=2**14)
    check.reorder('t2',first=False)
    forplot = abs(check)['t2':(-50,50)]
    myline = forplot.C.argmax('t2')
    forplot.smoosh(['ph1','ph2','power'],'indirect'
            ).setaxis('indirect','#'
                    ).reorder('indirect',first=False)
    if fl is not None:
        fl.push_marker()
        fl.next('cross-correlation')
        fl.image(forplot)
        fl.plot(myline.C.smoosh(['ph1','ph2','power'],'indirect'
            ).setaxis('indirect','#'
                    ).reorder('indirect',first=False)
                , 'w', linewidth=3, alpha=0.5)
    return s*exp(-1j*2*pi*myline*s.fromaxis('t2'))
