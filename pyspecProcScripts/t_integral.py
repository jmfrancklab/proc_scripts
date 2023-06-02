from pyspecdata import *
from pyspecProcScripts import *
import matplotlib.pyplot as plt
import numpy as np

def t_integral_w_error(
        convolved,
        original,
        sig_pathway,
        apo_fn,
        cutoff = 0.1,
        frq_slice = None,
        direct = 't2',
        fl = None,
        return_frq_slice = False):
    """Calculates the propagation of error for a convolved signal in the time domain (we know the data points in the time domain remain uncorrelated to each other). Returns the integrated signal of the convolved function with the correct error associated.

    You will need to feed both the convolved data as well as the data prior to convolution for this to correctly calculate the error.

    Parameters
    ==========
    convolved:      nddata
                    The data that has been apodized in the time domain and 
                    therefore convolved in the frequency domain
    original:       nddata
                    The data prior to convolution.
    sig_pathway:    dict
                    Dictionary of the path of the desired signal.
    apo_fn:         nddata
                    The function that you are apodizing with
    cutoff:         float
                    Value fed to the integrate_limits to decide how to cut the signal. 
                    The larger this value the narrower the cut as it decides the cutoff 
                    to be cutoff * max of the signal.
    frq_slice:      tuple
                    The user can feed a manual frq slice if they don't want to use the 
                    automated version or need to maintain a slice.
    direct:         str
                    Direct axis.
    
    Returns
    =======
    s:          nddata
                Integrated data with error associated with coherence pathways not included 
                in the signal pathway. variance is calculated in the time domain.
    frq_slice:  tuple
                the frequency slice that was automated using integrate_limits
    """            
    # Find signal slice and slice signal out
    assert convolved.get_ft_prop(direct), "need to be in frequency domain to decide the slice"
    assert original.get_ft_prop(direct), "need to be in frequency domain to decide the slice"
    if frq_slice is not None:
        frq_slice = frq_slice
    else:
        frq_slice = integrate_limits(
                select_pathway(convolved, sig_path),  
                cutoff = cutoff, fl=fl)
    logging.debug(strm("frq_slice is", frq_slice))  
    d = original.C
    s = convolved[direct:frq_slice].C
    #}}}
    #{{{calculate error
    apo_fn['t2':0] *= 0.5
    integral_apo = apo_fn.C.run(lambda x: x**2).integrate('t2')
    if fl is not None:
        fl.next('Time domain')
        fl.plot(select_pathway(d,signal_pathway,label = 'unapodized data',alpha = 0.5))
        fl.plot(apo_fn*abs(d.C).max(), label = 'apo function', alpha = 0.5)
    dt = np.diff(original.real.C.getaxis(direct)[r_[0,1]]).item()    
    temp = select_pathway(d,sig_pathway)
    temp.data[:] = nan
    if fl is not None:
        fl.next('masked out signal pathway')
        fl.image(d)
    t_var = np.var(d.real.data[np.isfinite(d.data)],ddof=1)
    var_t = (t_var * dt * integral_apo)
    std_t = sqrt(var_t.data)
    #}}}
    convolved.ift('t2')
    t_int = select_pathway(convolved.C,sig_pathway)
    retval = t_int.integrate(direct).set_error(np.asarray(std_t))
    if not return_frq_slice:
        return retval
    elif return_frq_slice:
        return retval, frq_slice
