import pyspecdata as psp
from .integrate_limits import integrate_limits
from .simple_functions import select_pathway
import logging
import numpy as np


def integral_w_errors(
    s,
    sig_path,
    error_path,
    cutoff = 0.25,
    indirect="vd",
    direct="t2",
    fl=None,
    return_frq_slice=False,
):
    """Calculates the propagation of error for the given signal and returns
    signal with the error associated.

    Before declaring the error_path,
    look at an examples such as integration_w_error.py to see how to decide which
    excluded pathways to take the error over.

    Parameters
    ==========
    sig_path:   dict
                Dictionary of the path of the desired signal.
    error_path: dict
                Dictionary of all coherence pathways that are
                not the signal pathway.

                Before declaring the error_path,
                look at an examples such as integration_w_error.py to see how to decide which
                excluded pathways to take the error over.
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
    frq_slice = integrate_limits(
        select_pathway(s, sig_path),  
        cutoff=cutoff, fl=fl)
    logging.debug(psp.strm("frq_slice is", frq_slice))
    s = s[direct:frq_slice]
    f = s.getaxis(direct)
    df = f[1] - f[0]
    errors = []
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect, direct])
    extra_dims = [j for j in all_labels if not j.startswith("ph")]
    if len(extra_dims) > 0:
        raise ValueError(
            "You have extra (non-phase cycling, non-indirect) dimensions: "
            + str(extra_dims))
    collected_variance = psp.ndshape([psp.ndshape(s)[indirect], len(error_path)], 
        [indirect, "pathways"]).alloc()
    for j in range(len(error_path)):
        # calculate N₂ Δf² σ², which is the variance of the integral (by error propagation)
        # where N₂ is the number of points in the indirect dimension
        s_forerror = select_pathway(s, error_path[j])
        # previous line wipes everything out and starts over -- why not use
        # collected_variance above, as I had originally set up --> part of
        # issue #44
        if j == 0:
            N2 = psp.ndshape(s_forerror)[direct]
        # mean divides by N₁ (indirect), integrate multiplies by Δf, and the
        # mean sums all elements (there are N₁N₂ elements)
        s_forerror -= s_forerror.C.mean_all_but([indirect])
        s_forerror.run(lambda x: abs(x) ** 2 / 2).mean_all_but([indirect])
        s_forerror *= df ** 2  # Δf
        s_forerror *= N2
    s = select_pathway(s, sig_path)
    retval = s.integrate(direct).set_error(psp.sqrt(s_forerror.data))
    if not return_frq_slice:
        return retval
    elif return_frq_slice:
        return retval, frq_slice


def active_propagation(
    s, signal_path, cutoff = 0.25, indirect="vd", direct="t2", fl=None, offset=500.0
):
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
    frq_slice = integrate_limits(
            select_pathway(s, signal_path),
            cutoff = cutoff, fl=fl)
    logging.debug(psp.strm("frq_slice is", frq_slice))
    slice_span = frq_slice[-1]-frq_slice[0]
    slice_start = frq_slice[-1] + offset 
    off_res_frq_slice = (slice_start,slice_start+slice_span)
    full_s = s.C
    s = s[direct : off_res_frq_slice]  # grab all data more than
    #                                             offset to the right of the
    #                                             peak
    f = s.getaxis(direct)
    df = f[1] - f[0]
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect, direct])
    extra_dims = [j for j in all_labels if not j.startswith("ph")]
    if len(extra_dims) > 0:
        raise ValueError(
            "You have extra (non-phase cycling, non-indirect) dimensions: "
            + str(extra_dims))
    s_forerror = select_pathway(s, signal_path)
    N2 = psp.ndshape(s_forerror)[direct]
    s_forerror -= s_forerror.C.mean_all_but([indirect])
    s_forerror.run(lambda x: np.real(x) ** 2).mean_all_but([indirect])
    s_forerror *= df ** 2
    s_forerror *= N2
    s = select_pathway(s,signal_path)
    retval = s.integrate(direct).set_error(psp.sqrt(s_forerror.data))
    return retval

