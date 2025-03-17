import pyspecdata as psd
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


def _masked_var_multi(x, var_has_imag=False, axis=None):
    """Calculates the variance along a 1D axis.
    If the data is complex you must assign var_has_imag as true.
    By default it calculates the variance of the real of the data"""
    assert axis is not None

    def masked_var(x):
        if var_has_imag:  # take average of variance along real and imag
            return np.var(x[np.isfinite(x)], ddov=1) / 2
        else:
            return np.var(x[np.isfinite(x)], ddof=1)

    return np.apply_along_axis(masked_var, axis, x)


# }}}
def calc_masked_error(
    s,
    frq_slice,
    signal_pathway,
    excluded_pathways = None,
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
        Frequency range over which the signal resides.
    signal_pathway:   dict
                Dictionary of the path of the desired signal.
    excluded_pathways: list
                List of tuples containing the coherence pathways that are
                to be masked out when calculating the error.
                If no excluded_pathways are fed, the function will apply a mask
                over the signal pathway as well as over the integraion slice in
                all other coherence transfer pathways.
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
    df = s.get_ft_prop(direct, "df")
    N = len(s[direct])
    #if excluded_pathways is not None:
    # {{{ Define the pathways used for calculating the error
    collected_variance = s.C  # so we don't alter s
    phcycdims = [j for j in s.dimlabels if j.startswith("ph")]
    temp = select_pathway(collected_variance, signal_pathway)
    temp.data[:] = np.nan
    if excluded_pathways is not None:
        error_paths = [{phcycdims[i]: path[i] for i in range(len(phcycdims))} for
                path in excluded_pathways]
        for j in range(len(error_paths)):
            temp = select_pathway(collected_variance,error_paths[j])
            temp.data[:] = np.nan
    else:
        temp = collected_variance[direct:frq_slice]
        temp.data[:] = np.nan
    if fl is not None:
        fl.next("Frequency Noise")
        fl.image(collected_variance)
    collected_variance.run(_masked_var_multi, direct)
    for j in [k for k in s.dimlabels if k.startswith("ph")]:
        collected_variance.run(_masked_mean_multi, j)
    collected_variance = collected_variance * df**2 * N
    return collected_variance
