"""
Time-Domain Noise
=================

Here, we want to calculate the time-domain variance to use in error
propagation.  But, to make sure we calculate only noise, we want to mask out
portions of the frequency domain.
We propose that if we use a unitary transform,
Parseval's theorem tells us we can calculate :math:`\sigma_t^2` (time domain
variance) directly from the frequency-domain variance, *i.e.*
:math:`\sigma_\nu^2=\sigma_t^2`.
To confirm this, we construct a "spectrum" of pure noise and generate a
frequency-masked noise, and show that is the same as the unmasked time-domain.
"""

import numpy as np
from numpy import r_
from pyspecdata import *
from pyspecProcScripts import *

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
direct_t_dom_std = sqrt(example_data.C.run(np.var, "t").mean("ph") / 2)
# the way that we do FT is parseval preserved?
temp = example_data.C.ft("t", shift=True)
dt = temp.get_ft_prop("t", "dt")
df = temp.get_ft_prop("t", "df")
sigma_nu = sqrt(temp.run(np.var, "t").mean("ph") / 2)
print(
    "If we apply just FT as we normally would the std in the frequency"
    " domain is:",
    sigma_nu,
    "\nNote that this does **not** match σₜ, which is 1\n",
    "\nWhile I could use a unitary FT, I don't usually do that, so instead, I"
    " scale to move from σ_ν to σₜ.  Note I multiply by the √Δν and divide by"
    " √Δt to move *to* σₜ",
    sigma_nu * sqrt(df / dt),
)

# now, I can just calculate sigma_nu, since it's easier to mask out regions of
# the coherence domain where I expect there is signal (or phase cycling noise).
# I can then convert that to sigma_t
example_data.ft("t",shift=True,unitary=True)
example_data.ft("ph")

# {{{ I'm doing a mildly odd thing where I'm using "nan" to identify signal I
#     want to exclude from the variance calculation -- i.e. to mask it.  This
#     is assuming that I have signal that I'm not interested in including in
#     the calculation.
temp = select_pathway(example_data, signal_pathway)
temp.data[:] = nan  # note how I am NOT acting on a copy -- I am trying to
#                    manipulate the data at its original memory position!
# for the most complicated case I'll also say I want to exclude phase cycling
# noise -- so also exclude everything from the signal bandwidth this will give
# a conservative (small) estimate of the noise
temp = example_data["t":signal_window]
temp.data[:] = nan
# }}}
with figlist_var(black=True) as fl:
    fl.next("show the mask in white")
    forplot = example_data.C
    # in pyspecdata, nan shows up as the opposite (black vs. white) color vs. 0
    fl.image(forplot)
# {{{ Calculate the variance using new functions
#    now, I can do this:
result = calc_masked_variance(example_data)
# convert from σ_ν to σ_t
result *= result.get_ft_prop("t", "df") / result.get_ft_prop("t", "dt")
result.run(
    lambda x: sqrt(x)
)  # convert variance to std for subsequent comparison
print("The std when using the mask on unitary data is:", example_data)
print(
    "Because we can use the mask in the DCCT domain to exclude signal, that is"
    " the number we will want, in general."
)
print(
    "However, here we know that all of our data is noise, and so we should"
    " make sure that this matches the naive, direct time-domain calculation."
)
print("If it does, all the following numbers will be about 1.0:")
print(result / direct_t_dom_std)
# }}}
