"""
Check integral error calculation
================================

Generate a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise.
Check that the following match:

- integral w/ error (the canned routine
  :func:`~pyspecProcScripts.integral_w_errors`)
- propagate error based off the programmed σ of the normal distribution
- set the error bars based on the standard deviation (along the repeats
  dimension) of the *real* part of the integral
- propagate error based off the variance of the noise in the inactive
  coherence channels (do this manually inside this script -- should mimic
  what :func:`~pyspecProcScripts.integral_w_errors` does)
"""

from numpy import diff, r_, sqrt, real, exp, pi
from pyspecdata import ndshape, nddata, init_logging, figlist_var
from pyspecProcScripts import integral_w_errors

# sphinx_gallery_thumbnail_number = 1

init_logging(level="debug")
fl = figlist_var()
t2 = nddata(r_[0:1:1024j], "t2")
vd = nddata(r_[0:1:40j], "vd")
ph1 = nddata(r_[0, 2] / 4.0, "ph1")
ph2 = nddata(r_[0:4] / 4.0, "ph2")
signal_pathway = {"ph1": 0, "ph2": 1}
excluded_pathways = [
    signal_pathway,
    {"ph1": 0, "ph2": 0},
    {"ph1": 0, "ph2": 3},
]
# this generates fake clean_data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 10 Hz
clean_data = (
    21 * (1 - 2 * exp(-vd / 0.2)) * exp(+1j * 2 * pi * 100 * t2 - t2 * 10 * pi)
)
clean_data *= exp(signal_pathway["ph1"] * 1j * 2 * pi * ph1)
clean_data *= exp(signal_pathway["ph2"] * 1j * 2 * pi * ph2)
clean_data["t2":0] *= 0.5
fake_data_noise_std = 2.0
clean_data.reorder(["ph1", "ph2", "vd"])
bounds = (0, 200)  # seem reasonable to me
result = 0
n_repeats = 100
all_results = ndshape(clean_data) + (n_repeats, "repeats")
all_results.pop("t2").pop("ph1").pop("ph2")
all_results = all_results.alloc()
all_results.setaxis("vd", clean_data.getaxis("vd"))
print("shape of all results", ndshape(all_results))
for j in range(n_repeats):
    data = clean_data.C
    data.add_noise(fake_data_noise_std)
    # at this point, the fake data has been generated
    data.ft(["ph1", "ph2"])
    # {{{ usually, we don't use a unitary FT -- this makes it unitary
    data /= 0.5 * 0.25  # the dt in the integral for both dims
    data /= sqrt(ndshape(data)["ph1"] * ndshape(data)["ph2"])  # normalization
    # }}}
    dt = diff(data.getaxis("t2")[r_[0, 1]]).item()
    data.ft("t2", shift=True)
    # {{{
    data /= sqrt(ndshape(data)["t2"]) * dt
    s_int, frq_slice = integral_w_errors(
        data,
        sig_path=signal_pathway,
        excluded_pathways=excluded_pathways,
        excluded_frqs=[bounds],
        indirect="vd",
        fl=fl,
        return_frq_slice=True,
    )
    # }}}
    manual_bounds = data["ph1", 0]["ph2", 1]["t2":frq_slice]
    N = ndshape(manual_bounds)["t2"]
    df = diff(data.getaxis("t2")[r_[0, 1]]).item()
    manual_bounds.integrate("t2")
    # N terms that have variance given by fake_data_noise_std**2 each
    # multiplied by df
    all_results["repeats", j] = manual_bounds
    print("#%d" % j)
std_off_pathway = (
    data["ph1", 0]["ph2", 0]["t2":bounds]
    .C.run(lambda x: abs(x) ** 2 / 2)  # sqrt2 so variance is variance of real
    .mean_all_but(["t2", "vd"])
    .mean("t2")
    .run(sqrt)
)
print(
    "off-pathway std", std_off_pathway, "programmed std", fake_data_noise_std
)
propagated_variance_from_inactive = N * df**2 * std_off_pathway**2
# removed factor of 2 in following, which shouldn't have been there
propagated_variance = N * df**2 * fake_data_noise_std**2
fl.next("different types of error")
fl.plot(s_int, ".", capsize=6, label="std from int w err", alpha=0.5)
manual_bounds.set_error(sqrt(propagated_variance))
fl.plot(
    manual_bounds,
    ".",
    capsize=6,
    label=r"propagated from programmed variance",
    alpha=0.5,
)
all_results.run(real).mean("repeats", std=True)
# by itself, that would give error bars, but the data would be
# averaged -- better to put the data in the same position
manual_bounds.set_error(all_results.get_error())
# the fact that this matches the previous shows that my sample size is
# large enough to give good statistics
fl.plot(manual_bounds, ".", capsize=6, label=r"std from repeats", alpha=0.5)
manual_bounds.set_error(sqrt(propagated_variance_from_inactive.data))
fl.plot(
    manual_bounds,
    ".",
    capsize=6,
    label=r"propagated from inactive std",
    alpha=0.5,
)
fl.show()
