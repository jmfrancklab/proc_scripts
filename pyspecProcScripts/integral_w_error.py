import pyspecdata as psp
from .integrate_limits import integrate_limits
from .simple_functions import select_pathway
from .calc_error import calc_masked_variance
import logging
import numpy as np


def frequency_domain_integral(
    s,
    signal_pathway=None,
    excluded_frqs=None,
    excluded_pathways=None,
    cutoff=0.1,
    convolve_method="Gaussian",
    indirect="vd",
    direct="t2",
    fl=None,
    return_frq_slice=False,
):
    """Returns the integral along the direct dimension,
    with automatically determined bounds.
    Error is calculated by propagation of error in the frequency domain,
    and included in the resulting nddata.

    Before declaring the error_path,
    look at an examples such as integration_w_error.py to see how to
    decide which excluded pathways to take the error over.

    Used in these examples:
    `examples/integration_with_error.py`
    `examples/check_integration_error.py`
    `examples/broken/test_error.py`

    Parameters
    ==========
    signal_pathway: dict
        Dictionary of the path of the desired signal.
    excluded_frqs: list
        List of tuples containing frequencies to be filtered out when
        calculating the variance of the spectral datapoints.

        If this is set to `None`, then the code will use the
        automatically chosen integration slice.

        If this is set to `[]` (empty list),
        then the code will not exclude any frequencies.
    excluded_pathways: list
        List of dictionaries containing all coherence pathways that are to be
        masked out when calculating the error. This should include the signal
        pathway!
    cutoff: float
        Multiplier used in determining the integral bounds. Higher values will
        produce smaller integrals.
    convolve_method: str
        method of convolution used in integrating limits passed on to
        :func:`integrate_limits`
    indirect: str
        Indirect axis.
    direct: str
        Direct axis.

    Returns
    =======
    s: nddata
        Data sampled at the indicated coherence pathway, and integrated
        over automatically determined bounds.
        The bounds are found by applying a matched filter and using the
        "cutoff" as a fraction of the maximum signal intensity to
        determine the limits.
        Errors are determined by error propagation in the frequency
        domain, where the noise associated with the spectral datapoints
        determined from standard deviation of the DCCT data,
        where the masked frequencies and coherence pathways have been
        excluded.
    """
    signal_pathway = signal_pathway or s.get_prop("coherence_pathway")
    assert s.get_ft_prop(direct), "need to be in frequency domain!"
    if convolve_method is not None:
        kwargs = {"convolve_method": convolve_method}
    else:
        kwargs = {}
    frq_slice = integrate_limits(
        select_pathway(s, signal_pathway),
        cutoff=cutoff,
        fl=fl,
        **kwargs,
    )
    df = s.get_ft_prop(direct, "df")
    # {{{ round to the nearest discrete coord
    for j in range(2):
        idx = np.searchsorted(s[direct], frq_slice[j] + 0.5 * df)
        frq_slice[j] = s[direct][idx]
    # }}}
    logging.debug(psp.strm("frq_slice is", frq_slice))
    if excluded_frqs is None:
        excluded_frqs = [frq_slice]
    spectral_datapoint_variance = calc_masked_variance(
        s,
        excluded_frqs=excluded_frqs,
        indirect=indirect,
        excluded_pathways=excluded_pathways,
        fl=fl,
    )
    s = s[direct:frq_slice]
    N = s.shape[direct]  # number of pts within the slice
    all_labels = set(s.dimlabels)
    all_labels -= set([indirect, direct])
    extra_dims = [j for j in all_labels if not j.startswith("ph")]
    if len(extra_dims) > 0:
        raise ValueError(
            "You have extra (non-phase cycling, non-indirect) dimensions: "
            + str(extra_dims)
        )
    # returns ∫s(ν)dν with error set to √(σ²_ν Δν² N)
    retval = (
        select_pathway(s, signal_pathway)
        .integrate(direct)
        .set_error(psp.sqrt(spectral_datapoint_variance.data * df**2 * N))
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
