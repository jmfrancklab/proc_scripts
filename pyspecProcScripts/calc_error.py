import pyspecdata as psd
import numpy as np
from .simple_functions import select_pathway


def calc_error(
    s,
    frq_slice,
    signal_pathway,
    error_path,
    direct="t2",
    indirect="nScans",
    use_mask=False,
    fl = None,
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
    df = s.get_ft_prop(direct,"df")
    collected_variance = (
        psd.ndshape(
            [psd.ndshape(s)[indirect], len(error_path)], [indirect, "pathways"]
        )
        .alloc()
        .setaxis("pathways", error_path)
    )
    if use_mask == False:
        assert error_path is not None, (
            "If you don't want to use the mask you need to tell me what"
            " pathways to use for calculating the noise!"
        )
        for j in range(len(error_path)):
            # calculate N₂ Δf² σ², which is the variance of the integral
            # (by error propagation) where N₂ is the number of points in
            # the indirect dimension
            s_forerror = select_pathway(s, error_path[j])
            # previous line wipes everything out and starts over -- why not use
            # collected_variance above, as I had originally set up --> part of
            # issue #44
            if j == 0:
                N2 = psd.ndshape(s_forerror)[direct]
            # mean divides by N₁ (indirect), integrate multiplies by Δf, and the
            # mean sums all elements (there are N₁N₂ elements)
            s_forerror -= s_forerror.C.mean_all_but([indirect, direct]).mean(
                direct
            )
            s_forerror.run(lambda x: abs(x) ** 2 / 2).mean_all_but(
                [direct, indirect]
            ).mean(direct)
            s_forerror *= df**2  # Δf
            s_forerror *= N2
            collected_variance["pathways", j] = s_forerror
        collected_variance = psd.sqrt(
            collected_variance.sum("pathways") / len(error_path)
        )
        return collected_variance
