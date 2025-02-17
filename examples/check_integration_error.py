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
- propagate error based off the variance of the noise in the inactive coherence
  channels (do this manually inside this script -- should mimic what
  :func:`~pyspecProcScripts.integral_w_errors` does)
"""
import pyspecdata as psd
import numpy as np
from numpy import r_
from pyspecProcScripts import integral_w_errors, select_pathway
from pylab import rcParams
import sympy as s
from collections import OrderedDict

rcParams["image.aspect"] = "auto"  # for sphinx gallery
# sphinx_gallery_thumbnail_number = 1

t2, vd, ph1, ph2 = s.symbols("t2 vd ph1 ph2")
n_repeats = 100
noise_bounds = (0, 200)  # seem reasonable to me
signal_pathway = {"ph1": 0, "ph2": 1}
# {{{ Generate fake IR dataset with multiple repeats
# This generates fake clean_data w/ a T₂ of 0.2s amplitude of 21, just to pick
# a random amplitude offset of 300 Hz, FWHM 10 Hz
expression = (
    21
    * (1 - 2 * s.exp(-vd / 0.2))
    * s.exp(+1j * 2 * s.pi * 100 * t2 - t2 * 10 * s.pi)
)
orderedDict = [
    ("vd", psd.nddata(r_[0:1:40j], "vd")),
    ("ph1", psd.nddata(r_[0, 2] / 4.0, "ph1")),
    ("ph2", psd.nddata(r_[0:4] / 4.0, "ph2")),
    ("t2", psd.nddata(r_[0:1:1024j], "t2")),
]
clean_data = psd.fake_data(
    expression, OrderedDict(orderedDict), signal_pathway, scale=0
)
clean_data.reorder(["vd", "t2"], first=False)
clean_data["t2":0] *= 0.5
fake_data_noise_std = 2.0
clean_data.reorder(["ph1", "ph2", "vd"])
# Unitary FT concentrates signal into a single index/position along coherence
# transfer pathway while preserve the vector norm so divide by sqrt(N_ph)
clean_data /= np.sqrt(
    psd.ndshape(clean_data)["ph1"] * psd.ndshape(clean_data)["ph2"]
)
# }}}
# Allocate nddata to place the calculated integrals with error into
all_results = psd.ndshape(clean_data) + (n_repeats, "repeats")
all_results.pop("t2").pop("ph1").pop("ph2")
all_results = all_results.alloc().setaxis("vd", clean_data.getaxis("vd"))
with psd.figlist_var() as fl:
    for j in range(n_repeats):
        data = clean_data.C  # need a copy so that each repeat is independent
        #                      and has it's own noise
        data.add_noise(fake_data_noise_std)  # add noise to each individual
        #                                      repeat
        # at this point, the fake data has been generated
        data.ft("t2")
        dt = data.get_ft_prop("t2", "dt")
        # {{{ Vector-normalize the FT
        data /= np.sqrt(psd.ndshape(data)["t2"]) * dt
        # }}}
        if j == 0:
            s_int, frq_slice = integral_w_errors(
                data,
                signal_pathway,
                fl=fl,
                return_frq_slice=True,
            )
        else:
            s_int, frq_slice = integral_w_errors(
                data,
                signal_pathway,
                fl=None,
                return_frq_slice=True,
            )
        manual_bounds = select_pathway(data["t2":frq_slice], signal_pathway)
        N = psd.ndshape(manual_bounds)["t2"]
        df = manual_bounds.get_ft_prop("t2", "df")
        manual_bounds.integrate("t2")
        # N terms that have variance given by fake_data_noise_std**2 each
        # multiplied by df
        all_results["repeats", j] = manual_bounds
        print("#%d" % j)
    std_off_pathway = (
        data["ph1", 0]["ph2", 0]["t2":noise_bounds]
        .C.run(
            lambda x: abs(x) ** 2 / 2
        )  # sqrt2 so variance is variance of real
        .mean_all_but(["t2", "vd"])
        .mean("t2")
        .run(np.sqrt)
    )
    print(
        "average of the off-pathway std: ",
        np.mean(std_off_pathway.data),
        "programmed std: ",
        fake_data_noise_std,
    )
    propagated_variance_from_inactive = N * df**2 * std_off_pathway**2
    propagated_variance = N * df**2 * fake_data_noise_std**2
    fl.next("different types of error")
    fl.plot(s_int, ".", capsize=6, label="std from int w err", alpha=0.5)
    manual_bounds.set_error(np.sqrt(propagated_variance))
    fl.plot(
        manual_bounds,
        ".",
        capsize=6,
        label=r"propagated from programmed variance",
        alpha=0.5,
    )
    all_results.run(np.real).mean("repeats", std=True)
    # by itself, that would give error bars, but the data would be averaged --
    # better to put the data in the same position
    manual_bounds.set_error(all_results.get_error())
    # the fact that this matches the previous shows that my sample size is
    # large enough to give good statistics
    fl.plot(
        manual_bounds, ".", capsize=6, label=r"std from repeats", alpha=0.5
    )
    manual_bounds.set_error(np.sqrt(propagated_variance_from_inactive.data))
    fl.plot(
        manual_bounds,
        ".",
        capsize=6,
        label=r"propagated from inactive std",
        alpha=0.5,
    )
