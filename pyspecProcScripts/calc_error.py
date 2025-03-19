from .simple_functions import select_pathway
import numpy as np


# {{{ Functions to calculate var and mean of nan-masked data
#     Pulled from time_domain_noise.py example
def _masked_mean_multi(x, axis=None):
    "Calculate the mean of nan-masked data on a 1D axis"
    assert axis is not None

    def masked_mean(x):
        return np.mean(x[np.isfinite(x)])

    return np.apply_along_axis(masked_mean, axis, x)


def _masked_var_multi(x, axis=None):
    """Calculates the variance along a 1D axis if an axis is fed or returns a 
    single float for the variance. The data must real"""

    def masked_var(x):
        if np.iscomplex(axis) and sum(abs(np.imag(axis))) > 1e-7:
            raise ValueError(
                "Data claims to be complex but has a negligible imaginary part"
            )
        elif not np.iscomplex(axis):
            return np.var(x[np.isfinite(x)], ddof=1)
        else:
            raise ValueError(
                "You're giving me complex data but I only calculate the"
                " variance of real data!"
            )

    if axis is not None:
        return np.apply_along_axis(masked_var, axis, x)
    else:
        return masked_var(x)


# }}}
def calc_masked_error(
    s,
    excluded_frqs=None,
    excluded_pathways=None,
    direct="t2",
    indirect="vd",
    fl=None,
):
    """Calculates the propagation of error for the given signal.

    Before declaring the excluded_pathways or excluded_frqs,
    look at an examples such as integration_w_error.py to see how to
    decide which excluded pathways to mask out.

    Parameters
    ==========
    s: nddata
        Full nddata set in the frequency domain for which the error is being
        propagated.
    excluded_frqs: list
        list of frequency ranges that will be filtered out in calculating the
        error - it is assumed this region in all coherence transfer pathways
        contains some amount of phase cycling noise.
    excluded_pathways: list
        List of dictionaries containing the coherence pathways that are
        to be masked out when calculating the error.
        If no excluded_pathways are fed, the function will apply a mask
        over just the signal pathway.
    direct: str
        Direct axis.
    indirect: str
        Indirect axis.
    Returns
    =======
    collected_variance: nddata
        The variance of the spectral datapoints with the only dimension being
        the indirect.
    """
    collected_variance = s.C  # so we don't alter s when we apply the mask
    # {{{ filter out excluded error pathways
    if isinstance(excluded_pathways, dict):
        raise ValueError(
            "excluded_pathways should be a list of dicts."
            "If you really mean to exclude only one pathway, pass a "
            "list with a single dict inside"
        )
    if excluded_pathways is None:
        try:
            excluded_pathways = [s.get_prop("coherence_pathway")]
        except Exception:
            assert excluded_pathways is not None, (
                "You need to at least give"
                " me the signal pathway to mask out, and your data does "
                "not have a property called coherence_pathway!"
            )
    if excluded_pathways is not None and len(excluded_pathways) > 0:
        for j in range(len(excluded_pathways)):
            temp = select_pathway(collected_variance, excluded_pathways[j])
            temp.data[:] = np.nan
    # }}}
    # {{{ Filter out frq_slice where ph noise resides
    if isinstance(excluded_frqs, tuple):
        raise ValueError(
            "excluded_frqs should be a list of tuples."
            "If you really mean to exclude only one slice, pass a "
            "list with a single tuple inside"
        )
    if excluded_frqs is not None and len(excluded_frqs) > 0:
        for j in range(len(excluded_frqs)):
            collected_variance[direct : excluded_frqs[j]].data[:] = np.nan
    # }}}
    if fl is not None:
        fl.next("Masked Frequency Noise")
        fl.image(collected_variance)
    # {{{ Average over remaining ct pathways so the returned nddata has one
    #     dimension
    for j in [k for k in s.dimlabels if k.startswith("ph")]:
        collected_variance.run(_masked_mean_multi, j)
    # }}}
    return collected_variance.run(_masked_var_multi, direct)
