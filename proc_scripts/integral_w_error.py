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
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect,direct])
    extra_dims = [j for j in all_labels if not j.startswith('ph')]
    if len(extra_dims) > 0:
        raise ValueError("You have extra (non-phase cycling, non-indirect) dimensions: "
                +str(extra_dims))
        print("LENGTH OF ERROR PATH:",len(error_path))
    collected_variance = ndshape(
         [ndshape(s)[indirect],len(error_path)],[indirect,'pathways']).alloc()
    check_var = collected_variance.C
    print("COLLECTED VARIANCE IS :",collected_variance)
    for j in range(len(error_path)):
        # calculate N₂ Δf² σ², which is the variance of the integral (by error propagation)
        # where N₂ is the number of points in the indirect dimension
        s_forerror = select_pathway(s,error_path[j])
        manual_bounds = s_forerror.C
        std_off_pathway = manual_bounds.C
        manual_bounds.integrate(direct)
        if j==0: N2 = ndshape(s_forerror)[direct]
        # mean divides by N₁ (indirect), integrate multiplies by Δf, and the
        # mean sums all elements (there are N₁N₂ elements)
        print("SHAPE OF S_FORERROR",ndshape(s_forerror))
        s_forerror.run(lambda x: abs(x)**2).mean_all_but([indirect,direct]).integrate(direct)
        s_forerror *= df # Δf
        check_var['pathways',j] = manual_bounds
        collected_variance['pathways',j] = s_forerror
        std_off_pathway = std_off_pathway.C.run(lambda x: abs(x)**2).mean_all_but([direct,indirect]).mean(direct).run(sqrt)
    collected_variance.mean('pathways') # mean the variance above across all pathways
    check_var.mean('pathways')
    # {{{ variance calculation for debug
    print("(inside automatic routine) the stdev seems to be",sqrt(collected_variance/(df*N2)))
    print("automatically calculated integral error:",sqrt(collected_variance.data))
    # }}}
    prop_var_from_inact = N2 *df**2*std_off_pathway**2
    s = select_pathway(s,sig_path)
    d = s.C
    if not return_frq_slice:
        return s.integrate(direct).set_error(sqrt(collected_variance.data))
    elif return_frq_slice:
        return s.integrate(direct).set_error(sqrt(collected_variance.data)), frq_slice, manual_bounds.set_error(sqrt(prop_var_from_inact.data))
