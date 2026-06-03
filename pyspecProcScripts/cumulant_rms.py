import numpy as np
from pyspecdata import mydiff


def cumulant_rms(d, indirect_dim, direct_dim="$B_0$"):
    """Return the normalized cumulative RMS change along an indirect axis.

    This diagnostic follows how much the ESR spectrum changes from one
    indirect point to the next.  Smooth alignment should reduce sharp
    horizontal-jitter artifacts, making this cumulant smoother along the
    indirect dimension.
    """
    # pyspecdata's .real property already returns a copy with a real-valued
    # view of the data, so an additional .C before .real is not needed.
    average_spectrum = d.real.C.mean(indirect_dim)
    cumulant = d.real.C.run(mydiff, indirect_dim)
    # mydiff preserves the axis length, but the final difference point is not
    # meaningful for this cumulant, so explicitly zero it before squaring.
    final_difference_point = [slice(None, None, None)]
    final_difference_point *= len(cumulant.dimlabels)
    final_difference_point[cumulant.dimlabels.index(indirect_dim)] = -1
    cumulant.data[tuple(final_difference_point)] = 0

    # Normalize the point-to-point differences by the vector norm of each
    # indirect trace, so the diagnostic emphasizes shape changes rather than
    # absolute signal level.
    cumulant /= d.C.run(lambda x: abs(x) ** 2).sum(direct_dim).run(np.sqrt)
    cumulant.run(lambda x: abs(x) ** 2).mean(direct_dim).run(np.sqrt).run(
        np.cumsum, indirect_dim
    )
    cumulant /= (
        average_spectrum.C.run(lambda x: abs(x) ** 2)
        .sum(direct_dim)
        .run(np.sqrt)
    )
    return cumulant
