from pyspecdata import *
import numpy as np

def integrate_limits(s, axis='t2', convwidth=100):
    frq_slice = s.C.mean_all_but(axis).convolve(axis,convwidth).contiguous(lambda x: abs(x)>0.5*abs(x).data.max())[0]
    return s[axis,frq_slice]
