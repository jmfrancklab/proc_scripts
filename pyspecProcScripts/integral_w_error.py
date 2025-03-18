import pyspecdata as psp
from .integrate_limits import integrate_limits
from .simple_functions import select_pathway
from .calc_error import calc_masked_error
import logging
import numpy as np


def integral_w_errors(
    s,
    sig_path,
    excluded_pathways=None,
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
    look at an examples such as integration_w_error.py to see how to
    decide which excluded pathways to take the error over.

    Parameters
    ==========
    sig_path:   dict
                Dictionary of the path of the desired signal.
    excluded_pathways: list
                List of tuples containing all coherence pathways that are
                to be masked out when calculating the error.
                By default the function will mask out the frequency slice of
                all coherence pathways and the entire signal pathway.
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
             Integrated data with error associated with coherence pathways
             not included in the signal pathway.
    """
    assert s.get_ft_prop(direct), "need to be in frequency domain!"
    if convolve_method is not None:
        kwargs = {"convolve_method": convolve_method}
    else:
        kwargs = {}
    frq_slice = integrate_limits(
        select_pathway(s, sig_path),
        cutoff=cutoff,
        fl=fl,
        **kwargs,
    )
    logging.debug(psp.strm("frq_slice is", frq_slice))
    # We will be calculating the error over the signal slice
    s = s[direct:frq_slice]
    # {{{ variables in calculating error over slice
    N = len(s[direct])  # number of pts within the slice
    df = s.get_ft_prop(direct, "df")
    # }}}
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect, direct])
    extra_dims = [j for j in all_labels if not j.startswith("ph")]
    if len(extra_dims) > 0:
        raise ValueError(
            "You have extra (non-phase cycling, non-indirect) dimensions: "
            + str(extra_dims)
        )
    variance = calc_masked_error(
        s,
        frq_slice,
        sig_path,
        indirect=indirect,
        excluded_pathways=excluded_pathways,
        fl=fl,
    )
    std_error = psp.sqrt(variance * df**2 * N)
    retval = (
        select_pathway(s, sig_path).integrate(direct).set_error(std_error.data)
    )
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
