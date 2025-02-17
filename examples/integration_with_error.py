"""
Check Integration_w_error Function
==================================

Generates a fake dataset of an inversion recovery w/ normally distributed
random noise.
Then checks that automatically chosen integral bounds perform similar to or
better than what you would choose by hand.
"""
import pyspecdata as psd
from pyspecProcScripts import integral_w_errors, select_pathway
from numpy.random import seed
import numpy as np
from numpy import r_
from pylab import rcParams
from collections import OrderedDict
import sympy as s

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2

seed(2021)  # so the same random result is generated every time -- 2021 is
#             meaningless
t2, vd, ph1, ph2 = s.symbols("t2 vd ph1 ph2")
signal_pathway = {"ph1": 0, "ph2": 1}
manual_slice = (60, 140)  # manually chosen integration bounds
# {{{ Generate fake IR dataset
# This generates fake data w/ a T₂ of 0.2s amplitude of 21, just to pick a
# random amplitude offset of 100 Hz, FWHM 10 Hz
data_expression = (
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
data = psd.fake_data(
    data_expression, OrderedDict(orderedDict), signal_pathway, scale=0
)
data.reorder(["vd", "t2"], first=False)
data["t2":0] *= 0.5
fake_data_noise_std = 2.0
data.reorder(["ph1", "ph2", "vd"])
# Unitary FT concentrates signal into a single index/position along coherence
# transfer pathway while preserve the vector norm so divide by sqrt(N_ph)
data /= np.sqrt(psd.ndshape(data)["ph1"] * psd.ndshape(data)["ph2"])
data.add_noise(fake_data_noise_std)
# }}}
data.ft("t2")
dt = data.get_ft_prop("t2", "dt")
# {{{ Vector-normalize the FT
data /= np.sqrt(psd.ndshape(data)["t2"]) * dt
# }}}
with psd.figlist_var() as fl:
    # {{{ First, run the code that automatically chooses integration bounds,
    #     integrates and assigns error
    s_int, returned_frq_slice = integral_w_errors(
        data, signal_pathway, fl=fl, return_frq_slice=True
    )
    fl.next("compare manual vs. automatic", legend=True)
    fl.plot(s_int, ".", label="fully auto: real", capsize=6)
    fl.plot(s_int.imag, ".", label="fully auto: imaginary", capsize=6)
    # }}}
    print("check the std after FT", np.std(data["ph1", 0]["ph2", 0].data.real))
    fl.next("Compare Manual vs. Automatic")
    # Run a controlled comparison between manually chosen integration bounds
    # and compare against automatically generated
    # Leave this as a loop so user can experiment with different bounds
    # As noted in issue #44 , manually chosen bounds underperform
    for bounds, thislabel in [
        (manual_slice, "manual bounds"),
        (tuple(returned_frq_slice), "auto bounds"),
    ]:
        # {{{ Calculate error of off-coherence pathway manually
        manual_bounds = select_pathway(data["t2":bounds], signal_pathway)
        assert manual_bounds.get_ft_prop("t2")
        # Check that the noise in an off coherence pathway correctly matches
        # the preset noise_std from above
        std_from_off_pathway = (
            data["ph1", 0]["ph2", 0]["t2":bounds]
            .C.run(lambda x: abs(x) ** 2 / 2)  # sqrt2 so variance is variance
            #                                   of real
            .mean_all_but(["t2", "vd"])
            .mean("t2")
            .run(np.sqrt)
        )
        print(
            "here is the average of the std calculated from an off pathway",
            np.mean(std_from_off_pathway.data),
            "does it match",
            fake_data_noise_std,
            "?",
        )
        # }}}
        # {{{ Manually calculate variance of data with manually set bounds
        N = psd.ndshape(manual_bounds)["t2"]
        df = manual_bounds.get_ft_prop("t2", "df")
        manual_bounds.integrate("t2")
        propagated_variance = N * df**2 * fake_data_noise_std**2
        propagated_variance_from_inactive = (
            N * df**2 * std_from_off_pathway**2
        )
        manual_bounds.set_error(np.sqrt(propagated_variance))
        fl.plot(
            manual_bounds,
            ".",
            capsize=6,
            label="%s (programmed σ)\n$%4g\\rightarrow%4g$"
            % ((thislabel,) + bounds),
            alpha=0.5,
        )
        manual_bounds.set_error(
            np.sqrt(propagated_variance_from_inactive.data)
        )
        fl.plot(
            manual_bounds,
            ".",
            capsize=6,
            label="%s (inactive CT σ)\n$%4g\\rightarrow%4g$"
            % ((thislabel,) + bounds),
            alpha=0.5,
        )
        # }}}
