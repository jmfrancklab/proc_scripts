import numpy as np
from pyspecdata import *
from pyspecProcScripts import *
#we know how to write a masked mean or std only along 1 dimension, so 
#   use numpy apply_alon_axis to make it a function that works along 1
#   dimension of multidimensional data

def masked_mean_multi(x, axis=None):
    "Calculates the mean on a 1D axis"
    assert axis is not None
    def masked_mean(x):
        "this only works for 1D data"
        return np.mean(x[np.isfinite(x)])
    return np.apply_along_axis(masked_mean,axis,x)
def masked_var_multi(x,axis=None, var_has_imag = True):
    "calculates the variance along a 1D axis"
    assert axis is not None
    def masked_var(x):
        "this only works for 1D data"
        if var_has_imag: # take average of variance along real and image
            return np.var(x[np.isfinite(x)], ddof=1)/2
        else:
            return np.var(x[np.isfinite(x)], ddof=1)
    return np.apply_along_axis(masked_var,axis,x)

