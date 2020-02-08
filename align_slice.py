from pyspecdata import *
def align_and_slice(s, dimension='t2', convwidth=500, threshold=0.05, fl=None):
    r"""Align the peak frequencies, and automatically slice them -- passing fl assumes you are debugging and want diagnostic plots"""
    if fl is not None:
        fl.push_marker()
    if not s.get_ft_prop(dimension):
        raise ValueError("Data has not been FTd yet!")
    s.ift(dimension)
    R = pi*convwidth # Lorentzian FWHM formula -- R = 1/T2, FWHM = R/pi
    s_forwidth = s.C*exp(-R * abs(s.fromaxis(dimension)))
    s_forwidth.ft(dimension)
    temp = s_forwidth.C.mean(dimension)
    temp /= abs(temp) 
    temp = s_forwidth / temp # all zeroth-order phased
    temp /= temp.C.real.sum(dimension) # normalize
    temp.run(real)
    if fl is not None:
        fl.next('diagnose align_and_slice -- normed and phased')
        fl.image(temp)
    temp *= temp.fromaxis(dimension)
    center_freq = temp.real.sum(dimension)
    print("center frequencies:",center_freq)
    s *= exp(-1j*2*pi*center_freq*s.fromaxis(dimension))
    s.ft(dimension)
    s_forwidth.mean_all_but([dimension])
    slices = s_forwidth.contiguous(lambda x:
            abs(x)>threshold*abs(s_forwidth.data).max())
    if fl is None:
        assert slices.shape[0] == 1, "found more than one peak: "+str(slices)
    slices = tuple(slices[0,:])
    print("slicing to",slices)
    if fl is not None:
        fl.pop_marker()
    return s[dimension:slices]
