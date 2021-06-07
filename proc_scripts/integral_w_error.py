from pyspecdata import *
from .integrate_limits import integrate_limits
import numpy as np
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval    
def integral_w_errors(s,sig_path,error_path, indirect='vd', direct='t2',fl=None,return_frq_slice=False):
    """Calculates the propagation of error for the given signal and returns
    signal with the error associated.
    
    Parameters
    ==========
    sig_path:   dict
                dictionary of the path of the desired signal
    error_path: dict
                dictionary of all coherence pathways that are 
                not the signal pathway
    indirect:   str
                indirect axis
    direct:     str
                direct axis
    
    Returns
    =======
    s:       nddata
                data with error associated with coherence pathways
                not included in the signal pathway
    """
    frq_slice = integrate_limits(select_pathway(s,sig_path),fl=fl)
    logging.debug(strm('frq_slice is',frq_slice))
    s = s[direct:frq_slice]
    f = s.getaxis(direct)
    df = f[1]-f[0]
    errors = []
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect,direct])
    extra_dims = [j for j in all_labels if not j.startswith('ph')]
    if len(extra_dims) > 0:
     raise ValueError("You have extra (non-phase cycling, non-indirect) dimensions: "
             +str(extra_dims))
    collected_variance = ndshape(
         [ndshape(s)[indirect],len(error_path)],[indirect,'pathways']).alloc()
    for j in range(len(error_path)):
     # calculate N₂ Δf² σ², which is the variance of the integral (by error propagation)
     # where N₂ is the number of points in the indirect dimension
     s_forerror = select_pathway(s,error_path[j])
     if j==0: N2 = ndshape(s_forerror)[direct]
     # mean divides by N₁ (indirect), integrate multiplies by Δf, and the
     # mean sums all elements (there are N₁N₂ elements)
     s_forerror -= s_forerror.C.mean_all_but([indirect, direct]).mean(direct)
     s_forerror.run(lambda x: abs(x)**2).mean_all_but([indirect,direct]).mean(direct)
     s_forerror *= df**2 # Δf
     s_forerror *= N2
     collected_variance['pathways',j] = s_forerror
    collected_variance.mean('pathways',std=True) # mean the variance above across all pathways
    # {{{ variance calculation for debug
    #print("(inside automatic routine) the stdev seems to be",sqrt(collected_variance/(df*N2)))
    print("automatically calculated integral error:",sqrt(collected_variance.data))
    # }}}
    s = select_pathway(s,sig_path)
    if not return_frq_slice:
        return s.integrate(direct).set_error(sqrt(collected_variance.data))
    elif return_frq_slice:
        return s.integrate(direct).set_error(sqrt(collected_variance.get_error())), frq_slice, sqrt(collected_variance.get_error())
