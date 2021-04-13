from pyspecdata import *
from .integrate_limits import integrate_limits
import numpy as np
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval    
def integral_w_errors(self,sig_path,error_path, indirect='vd', direct='t2'):
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
    self:       nddata
                data with error associated with coherence pathways
                not included in the signal pathway
    """
    frq_slice = integrate_limits(self)
    s = self[direct:frq_slice]
    t = self.getaxis(direct)
    dt = t[1]-t[0]
    errors = []
    all_labels = set(self.dimlabels)
    all_labels -= set([indirect,direct])
    extra_dims = [j for j in all_labels if not j.startswith('ph')]
    if len(extra_dims) > 0:
     raise ValueError("You have extra (non-phase cycling, non-indirect) dimensions: "
             +str(extra_dims))
    collected_variance = ndshape(
         [ndshape(s)['vd'],len(error_path)],['vd','pathways']).alloc()
    for j in range(len(error_path)):
     # calculate N₂ Δt² σ², which is the variance of the integral (by error propagation)
     # where N₂ is the number of points in the indirect dimension
     s_forerror = select_pathway(s,error_path[j])
     # mean divides by N₁ (indirect), integrate multiplies by Δt, and the
     # mean sums all elements (there are N₁N₂ elements)
     s_forerror.run(lambda x: abs(x)**2).mean_all_but([indirect,direct]).integrate(direct)
     s_forerror *= dt # Δt
     collected_variance['pathways',j] = s_forerror
    collected_variance.mean('pathways') # mean the variance above across all pathways
    s = select_pathway(s,sig_path)
    return s.integrate(direct).set_error(sqrt(collected_variance.data))
