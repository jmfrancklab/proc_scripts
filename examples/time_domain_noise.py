"""
Time-Domain Noise
=================

Here, we want to calculate the time-domain variance to use in error
propagation.
But, to make sure we calculate only noise, we want to mask out
portions of the frequency domain.
We propose that if we use a unitary transform,
Parseval's theorem
tells us we can calculate :math:`\\sigma_t^2` (time domain variance)
directly from the frequency-domain variance, *i.e.*
:math:`\\sigma_\\nu^2=\\sigma_t^2`.
To confirm this, we construct a "spectrum" of pure noise
and generate a frequency-masked noise,
and show that is the same as the unmasked time-domain noise.
"""

import numpy as np

from numpy import r_, sqrt, var
from pyspecdata import nddata, figlist_var
from pyspecProcScripts import calc_masked_variance

N = 1024
n_repeats = 50
signal_window = (-100, 200)  # wherever my "peak" shows up


# {{{ generate data with just noise with a phase cycling dimension and repeats
#     dimension
signal_pathway = {"ph": 1}
example_data = nddata(
    np.random.normal(size=4 * n_repeats * N)
    + 1j * np.random.normal(size=4 * n_repeats * N),
    [4, n_repeats, N],
    ["ph", "repeats", "t"],
)
example_data.setaxis("ph", r_[0:4] / 4)
example_data.setaxis("repeats", r_[0:n_repeats])
example_data.setaxis("t", r_[0 : 1 : 1j * N])
# }}}
# calculate the variance directly in the time domain.
# Because the data has no signal, know that this actually corresponds to the
# noise level:
direct_t_dom_std = (
    example_data.C.run(lambda x, axis=None: var(x, axis=axis, ddof=1) / 2, "t")
    .mean("ph")
    .mean("repeats")
    .run(sqrt)
)
# the way that we do FT is parseval preserved?
example_data.ft("t", shift=True)
dt = example_data.get_ft_prop("t", "dt")
df = example_data.get_ft_prop("t", "df")
sigma_nu = (
    example_data.C.run(lambda x, axis=None: var(x, axis=axis, ddof=1) / 2, "t")
    .mean("ph")
    .mean("repeats")
    .run(sqrt)
)
print(
    "If we apply just FT as we normally would the std in the frequency"
    " domain is:",
    sigma_nu,
    f"\nNote that this does **not** match σ_t, which is {direct_t_dom_std}\n",
    "\nWhile I could use a unitary FT, I don't usually do that, so instead, I"
    " scale to move from σ_ν to σ_t.  Note I multiply by the √Δν and divide by"
    " √Δt to move *to* σ_t, and I get",
    sigma_nu * sqrt(df / dt),
)

# now, I can just calculate sigma_nu,
# since it's easier to mask out regions of the coherence
# domain where I expect there is signal (or phase cycling noise).
# I can then convert that to sigma_t
example_data.ft("ph", unitary=True)

with figlist_var(black=True) as fl:
    # {{{ Calculate the variance using new functions
    #    now, I can do this:
    result = calc_masked_variance(
        example_data,
        direct="t",
        excluded_frqs=[signal_window],
        excluded_pathways=[signal_pathway],
        fl=fl,
    )
print(
    "The std when using the mask is:",
    sqrt(result.item()),
    "because there is no signal here, it should match σ_ν:",
    sigma_nu,
)
# next, convert from σ_ν to σ_t
result *= result.get_ft_prop("t", "df") / result.get_ft_prop("t", "dt")
result.run(sqrt)  # σ²_t→σ_t
print(
    "Because we can use the mask in the DCCT domain to exclude signal,"
    " we typically want to calculate σ_ν with a mask,"
    " and then convert to σ_t = σ_ν √(Δν/Δt)."
)
print(
    "Here, because we have all noise, we can check that our math works by"
    " seeing if the following number is about 1.0:"
)
print(result / direct_t_dom_std)
# }}}
