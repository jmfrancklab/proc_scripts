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
    """Calculates the variance along a 1D axis.
    If the data is complex you must assign var_has_imag as true
    so that the calculated variance is divided by 2.
    By default it calculates the variance of the real of the data"""

    def masked_var(x):
        if np.iscomplex(axis) and sum(abs(np.imag(axis))) > 1e-7:
            # take average of variance along real and imag
            return np.var(x[np.isfinite(x)], ddov=1) / 2
        else:
            return np.var(x[np.isfinite(x)], ddof=1)

    if axis is not None:
        return np.apply_along_axis(masked_var, axis, x)
    else:
        return masked_var(x)


# }}}
def calc_masked_error(
    s,
    frq_slice,
    signal_pathway,
    excluded_pathways=None,
    direct="t2",
    indirect="nScans",
    fl=None,
):
    """Calculates the propagation of error for the given signal.

    Before declaring the error_path,
    look at an examples such as integration_w_error.py to see how to
    decide which excluded pathways to take the error over.

    Parameters
    ==========
    s: nddata
        Full nddata set in the frequency domain for which the error is being
        propagated.
    frq_slice: tuple
        Frequency range that will be filtered out in calculating the error - it
        is assumed this region in all coherence transfer pathways contains some
        amount of phase cycling noise.
    excluded_pathways: list
                List of dictionaries containing the coherence pathways that are
                to be masked out when calculating the error.
                If no excluded_pathways are fed, the function will apply a mask
                over just the signal pathway.
    direct:     str
                Direct axis.
    indirect:   str
                Indirect axis.
    Returns
    =======
    collected_variance:       nddata
             The error associated with coherence pathways not included in the
             signal pathway.
    """
    collected_variance = s.C  # so we don't alter s
    # {{{ filter out signal pathway and excluded error pathways
    temp = select_pathway(collected_variance, signal_pathway)
    temp.data[:] = np.nan
    if type(excluded_pathways) == dict:
        raise ValueError(
            "excluded_pathways should be a list of dicts."
            "If you really mean to exclude only one pathway, pass a "
            "list with a single dict inside"
        )
    if excluded_pathways is not None and len(excluded_pathways) > 0:
        for j in range(len(excluded_pathways)):
            temp = select_pathway(collected_variance, excluded_pathways[j])
            temp.data[:] = np.nan
    # }}}
    # Filter out frq_slice where ph noise resides
    collected_variance[direct:frq_slice].data[:] = np.nan
    if fl is not None:
        fl.next("Frequency Noise")
        fl.image(collected_variance)
    for j in [k for k in s.dimlabels if k.startswith("ph")]:
        collected_variance.run(_masked_mean_multi, j)
    return _masked_var_multi(collected_variance.data).item()
