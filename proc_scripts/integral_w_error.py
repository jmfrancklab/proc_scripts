from pyspecdata import *
from .integrate_limits import integrate_limits
import numpy as np
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval    
def integral_w_errors(s,sig_path,error_path):
     frq_slice = integrate_limits(s)
     s_signal = s['t2':frq_slice].C
     s_signal = select_pathway(s,sig_path)
     errors = []
     for j in range(len(error_path)):
         s_forerror = s['t2':frq_slice]
         s_forerror = select_pathway(s_forerror,error_path[j])
         s_forerror.run(lambda x: abs(x)**2).mean_all_but(['vd','t2']).integrate('t2')
         errors.append(s_forerror)
     s_forerror = sum(errors)/len(errors)    
     sqrt(s_forerror)
     return s_signal.integrate('t2').set_error(s_forerror.data)
