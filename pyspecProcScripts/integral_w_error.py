import pyspecdata as psp
from .integrate_limits import integrate_limits
from .simple_functions import select_pathway
import logging
import numpy as np
def integral_w_errors(s,sig_path,error_path, convolve_method=None, indirect='vd', direct='t2',fl=None,return_frq_slice=False):
    """Calculates the propagation of error for the given signal and returns
    signal with the error associated.
    
    Parameters
    ==========
    sig_path:   dict
                Dictionary of the path of the desired signal.
    error_path: dict
                Dictionary of all coherence pathways that are 
                not the signal pathway.
    convolve_method: str
                method of convolution used in integrating limits
                passed on to :func:`integrate_limits`
    indirect:   str
                Indirect axis.
    direct:     str
                Direct axis.
    
    Returns
    =======
    s:       nddata
             Data with error associated with coherence pathways
             not included in the signal pathway.
    """
    assert s.get_ft_prop(direct), "need to be in frequency domain!"
    if convolve_method is not None:
        kwargs = {'convolve_method':convolve_method}
    else:
        kwargs = {}
    frq_slice = integrate_limits(select_pathway(s,sig_path),
            **kwargs)
    logging.debug(psp.strm('frq_slice is', frq_slice))
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
    collected_variance = psp.ndshape(
         [psp.ndshape(s)[indirect], len(error_path)],[indirect, 'pathways']).alloc()
    avg_error = []
    for j in range(len(error_path)):
        # calculate N₂ Δf² σ², which is the variance of the integral (by error propagation)
        # where N₂ is the number of points in the indirect dimension
        s_forerror = select_pathway(s,error_path[j])
        # previous line wipes everything out and starts over -- why not use
        # collected_variance above, as I had originally set up --> part of
        # issue #44 
        if j==0: N2 = psp.ndshape(s_forerror)[direct]
        # mean divides by N₁ (indirect), integrate multiplies by Δf, and the
        # mean sums all elements (there are N₁N₂ elements)
        s_forerror -= s_forerror.C.mean_all_but([indirect, direct]).mean(direct)
        s_forerror.run(lambda x: abs(x)**2/2).mean_all_but([direct,indirect]).mean(direct)
        s_forerror *= df**2 # Δf
        s_forerror *= N2
        avg_error.append(s_forerror)
    avg_error = sum(avg_error)/len(avg_error)
    # {{{ variance calculation for debug
    #print("(inside automatic routine) the stdev seems to be",sqrt(collected_variance/(df*N2)))
    #print("automatically calculated integral error:",sqrt(collected_variance.data))
    # }}}
    s = select_pathway(s,sig_path)
    retval = s.integrate(direct).set_error(psp.sqrt(s_forerror.data))
    if not return_frq_slice:
        return retval
    elif return_frq_slice:
        return retval, frq_slice

def active_propagation(s, signal_path, indirect='vd', direct='t2',fl=None,offset=500.0):
    """propagate error from the region `offset` to the right of the peak (where
    we assume there is only noise),  in the signal pathway `signal_path`, which
    we assume is the active coherence pathway.
    Include only the real part of the signal.

    Parameters
    ==========
    signal_path: dict
        Dictionary givin the active CT pathway
    indirect: str
        Name of the indirect dimension -- used to check that you don't have
        directions that are not direct, indirect, or phase cycling.
    direct: str
        Name of the direct dimension
    offset: float
        Distance (in Hz) between the auto-chosen integration bounds from
        :func:`integrate_limits` and the start of the "noise region."

    Returns
    =======
    retval: nddata
        just a data object with the error that this method predicts
    """
    assert s.get_ft_prop(direct), "need to be in frequency domain!"
    frq_slice = integrate_limits(select_pathway(s,signal_path),fl=fl)
    logging.debug(psp.strm('frq_slice is', frq_slice))
    s = s[direct:((frq_slice[-1]+offset),None)] # grab all data more than
    #                                             offset to the right of the
    #                                             peak
    df = s.get_ft_prop(direct,'df')
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect,direct])
    extra_dims = [j for j in all_labels if not j.startswith('ph')]
    if len(extra_dims) > 0:
        raise ValueError("You have extra (non-phase cycling, non-indirect) dimensions: "
                +str(extra_dims))
    s_forerror = select_pathway(s, signal_path)
    N = psp.ndshape(s_forerror)[direct]
    s_forerror.run(np.real).run(lambda x: abs(x)**2).mean_all_but([direct,indirect]).mean(direct)
    s_forerror *= df**2
    s_forerror *= N
    return s_forerror.run(psp.sqrt)
