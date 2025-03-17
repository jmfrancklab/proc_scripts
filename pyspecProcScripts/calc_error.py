import pyspecdata as psd
from .simple_functions import select_pathway
import numpy as np


# {{{ Functions to calculate var and mean of nan-masked data
#     Pulled from time_domain_noise.py example
def masked_mean_multi(x, axis=None):
    "Calculate the mean of nan-masked data on a 1D axis"
    assert axis is not None

    def masked_mean(x):
        return np.mean(x[np.isfinite(x)])

    return np.apply_along_axis(masked_mean, axis, x)


def masked_var_multi(x, var_has_imag=False, axis=None):
    "Calculates the variance along a 1D axis"
    assert axis is not None

    def masked_var(x):
        if var_has_imag:  # take average of variance along real and imag
            return np.var(x[np.isfinite(x)], ddov=1) / 2
        else:
            return np.var(x[np.isfinite(x)], ddof=1)

    return np.apply_along_axis(masked_var, axis, x)


# }}}
def calc_error(
    s,
    frq_slice,
    signal_pathway,
    error_path,
    direct="t2",
    indirect="nScans",
    use_mask=False,
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
    error_path: dict
                Dictionary of all coherence pathways that are
                not the signal pathway.

                Before declaring the error_path,
                look at an examples such as integration_w_error.py to
                see how to decide which excluded pathways to take the
                error over.
    direct:     str
                Direct axis.
    indirect:   str
                Indirect axis.
    use_mask: boolean
        If you would rather calculate the error excluding all the phase cycling
        noise this should be set to true. Specifically this will apply a mask
        over the signal pathway as well as over the integration slice in all
        other pathways. When False, the error will only be calculated for the
        pathways in error_path.

    Returns
    =======
    collected_variance:       nddata
             The error associated with coherence pathways not included in the
             signal pathway.
    """

    s = s[direct:frq_slice]
    df = s.get_ft_prop(direct, "df")
    N = len(s[direct])
    collected_variance = (
        psd.ndshape(
            [psd.ndshape(s)[indirect], len(error_path)], [indirect, "pathways"]
        )
        .alloc()
        .setaxis("pathways", error_path)
    )
    if use_mask is False:
        assert error_path is not None, (
            "If you don't want to use the mask you need to tell me what"
            " pathways to use for calculating the noise!"
        )
        for j in range(len(error_path)):
            # calculate N₂ Δf² σ², which is the variance of the integral
            # (by error propagation) where N₂ is the number of points in
            # the indirect dimension
            s_forerror = select_pathway(s, error_path[j])
            if j == 0:
                Nshape = psd.ndshape(s_forerror)[direct]
            # mean divides by N₁ (indirect), integrate multiplies by Δf, and
            # the mean sums all elements (there are N₁N₂ elements)
            s_forerror -= s_forerror.C.mean_all_but([indirect, direct]).mean(
                direct
            )
            s_forerror.run(lambda x: abs(x) ** 2 / 2).mean_all_but(
                [direct, indirect]
            ).mean(direct)
            s_forerror *= df**2  # Δf
            s_forerror *= Nshape
            collected_variance["pathways", j] = s_forerror
        collected_variance = psd.sqrt(
            collected_variance.sum("pathways") / len(error_path)
        )
    else:
        collected_variance = s.C  # so we don't alter s
        collected_variance[direct:frq_slice] = np.nan
        temp = select_pathway(collected_variance, signal_pathway)
        temp.data[:] = np.nan
        if fl is not None:
            fl.next("Frequency Noise")
            fl.image(collected_variance)
        collected_variance.run(masked_var_multi, direct)
        for j in [k for k in s.dimlabels if k.startswith("ph")]:
            collected_variance.run(masked_mean_multi, j)
        collected_variance = psd.sqrt(collected_variance * df**2 * N)
    return collected_variance
