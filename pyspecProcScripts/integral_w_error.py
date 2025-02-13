import pyspecdata as psp
from .integrate_limits import integrate_limits
from .simple_functions import select_pathway
import pyspecProcScripts as psdpr
import logging
import numpy as np
from numpy import r_

def masked_mean_multi(x, axis=None):
    "Calculates the mean on a 1D axis"
    assert axis is not None
    def masked_mean(x):
        "this only works for 1D data"
        return np.mean(x[np.isfinite(x)])
    return np.apply_along_axis(masked_mean,axis,x)
def masked_var_multi(x, var_has_imag=False, axis=None):
    "calculates the variance along a 1D axis"
    assert axis is not None
    def masked_var(x):
        "this only works for 1D data"
        if var_has_imag:
            return np.var(x[np.isfinite(x)], ddof=1) / 2
        else:
            return np.var(x[np.isfinite(x)],ddof=1)
    return np.apply_along_axis(masked_var, axis, x)
def integral_w_errors(
    s,
    sig_path,
    error_path,
    cutoff=0.1,
    convolve_method="Gaussian",
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
    if convolve_method is not None:
        kwargs = {"convolve_method": convolve_method}
    else:
        kwargs = {}
    frq_slice = integrate_limits(
        select_pathway(s, sig_path),
        convolve_method=convolve_method,
        cutoff=cutoff,
        fl=fl,
    )
    logging.debug(psp.strm("frq_slice is", frq_slice))
    noise = s.C
    noise[direct:frq_slice] = np.nan
    temp = psdpr.select_pathway(noise, sig_path)
    temp.data[:] = np.nan
    if fl is not None:
        fl.next("frequency noise masked")
        fl.image(noise)
    noise.run(masked_var_multi, direct)
    for j in [k for k in s.dimlabels if k.startswith("ph")]:
        noise.run(masked_mean_multi,j)
    s_int = psdpr.select_pathway(s,sig_path)[direct:frq_slice]    
    df = s_int.get_ft_prop(direct,"df")
    N = psp.ndshape(s)[direct]
    retval = s_int.integrate(direct).set_error(psp.sqrt(noise.data * df**2 * N))
    if not return_frq_slice:
        return retval
    elif return_frq_slice:
        return retval, frq_slice


def active_propagation(
    s, signal_path, indirect="vd", direct="t2", fl=None, offset=500.0
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
    frq_slice = integrate_limits(select_pathway(s, signal_path), fl=fl)
    logging.debug(psp.strm("frq_slice is", frq_slice))
    s = s[direct : ((frq_slice[-1] + offset), None)]  # grab all data more than
    #                                             offset to the right of the
    #                                             peak
    df = s.get_ft_prop(direct, "df")
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect, direct])
    extra_dims = [j for j in all_labels if not j.startswith("ph")]
    if len(extra_dims) > 0:
        raise ValueError(
            "You have extra (non-phase cycling, non-indirect) dimensions: "
            + str(extra_dims)
        )
    s_forerror = select_pathway(s, signal_path)
    N = psp.ndshape(s_forerror)[direct]
    s_forerror.run(np.real).run(lambda x: abs(x) ** 2).mean_all_but(
        [direct, indirect]
    ).mean(direct)
    s_forerror *= df**2
    s_forerror *= N
    return s_forerror.run(psp.sqrt)
