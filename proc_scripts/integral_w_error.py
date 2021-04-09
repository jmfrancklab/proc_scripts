from pyspecdata import *
from .integrate_limits import integrate_limits
import numpy as np
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval    
def integral_w_errors(self,sig_path,error_path, indirect='vd'):
     frq_slice = integrate_limits(self)
     s = self['t2':frq_slice]
     t = self.getaxis('t2')
     dt = t[1]-t[0]
     errors = []
     all_labels = set(self.dimlabels)
     all_labels -= set([indirect])
     extra_dims = [j for j in all_labels if not j.startswith('ph')]
     if len(extra_dims) > 0:
         raise ValueError("You have extra (non-phase cycling, non-indirect) dimensions: "
                 +str(extra_dims))
     collected_variance = ndshape(
             [ndshape(s)['vd'],len(error_path)],['vd','pathways']).alloc()
     for j in range(len(error_path)):
         s_forerror = select_pathway(s,error_path[j])
         s_forerror.run(lambda x: abs(x)**2).mean_all_but([indirect,'t2']).integrate('t2')
         s_forerror *= dt
         collected_variance['pathways',j] = s_forerror
     collected_variance.mean('pathways')
     s = select_pathway(s,sig_path)
     return s.integrate('t2').set_error(collected_variance.data)
