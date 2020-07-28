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
