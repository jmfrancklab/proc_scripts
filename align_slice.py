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
    s_conv.ft(dimension)
    first_time = True
    total_center_frq = 0
    j = 0
    while first_time or any(center_freq.data > 0.5):
        print("iterating")
        fl.next("iteration %d"%j)
        temp = abs(s_conv)
        fl.image(temp)
        temp /= temp.C.sum(dimension)
        temp *= temp.fromaxis(dimension)
        center_freq = temp.sum(dimension)
        print("center frequencies:",center_freq)
        s_conv.ift(dimension)
        s_conv *= exp(-1j*2*pi*center_freq*s.fromaxis(dimension))
        s_conv.ft(dimension)
        total_center_frq += center_freq
        j += 1
        if j > 10: break
    print("total frq:",total_center_frq)
    s *= exp(-1j*2*pi*total_center_frq*s.fromaxis(dimension))
    s.ft(dimension)
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
