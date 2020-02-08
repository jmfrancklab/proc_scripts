from pyspecdata import *
def align_and_slice(s, dimension='t2', convwidth=500):
    if not s.get_ft_prop(dimension):
        raise ValueError("Data has not been FTd yet!")
    temp = s.C.mean(dimension)
    temp /= abs(temp)
    temp *= s # all zeroth-order phased
    temp /= temp.C.sum(dimension) # normalize
    temp *= temp.fromaxis(dimension)
    center_freq = temp.real.sum(dimension)
    s.ift(dimension)
    s *= exp(-1j*2*pi*center_freq*s.fromaxis(dimension))
    R = 1./pi/convwidth # Lorentzian FWHM formula
    s_forwidth = s.C*exp(-R*s.getaxis(dimension))
    s.ft(dimension)
    s_forwidth.ft(dimension)
    ph0 = s_forwidth.C.mean(dimension)
    ph0 /= abs(ph0)
    s_forwidth /= ph0
    s_forwidth.mean_all_but([dimension])
    slices = s_forwidth.contiguous(lambda x:
            abs(x)>0.05*abs(s_forwidth.data).max())
    assert slices.shape[0] == 1, "found more than one peak"
    slices = tuple(slices[0,:])
    print("slicing to",slices)
    return s[dimension:slices]
