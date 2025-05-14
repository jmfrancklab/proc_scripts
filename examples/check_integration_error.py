"""
Check integral error calculation
================================

Generate a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise.
Check that the following match:

- integral w/ error (the canned routine
  :func:`~pyspecProcScripts.frequency_domain_integral`)
- propagate error based off the programmed σ of the normal distribution
- set the error bars based on the standard deviation (along the repeats
  dimension) of the *real* part of the integral
- propagate error based off the variance of the noise in the inactive
  coherence channels (do this manually inside this script -- should mimic
  what :func:`~pyspecProcScripts.frequency_domain_integral` does)
"""

# TODO ☐: read through current version, and make sure that the
#         description above matches

from numpy import diff, r_, sqrt, real, exp, pi, var
from pyspecdata import ndshape, nddata, init_logging, figlist_var
from pyspecProcScripts import frequency_domain_integral, select_pathway
from pylab import seed

# sphinx_gallery_thumbnail_number = 1
seed(2021)
init_logging(level="debug")
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
manual_bounds = (0, 200.0)
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
result = 0
n_repeats = 100
all_results_manual = ndshape(clean_data) + (n_repeats, "repeats")
all_results_manual.pop("t2").pop("ph1").pop("ph2")
all_results_manual = all_results_manual.alloc()
all_results_manual.setaxis("vd", clean_data.getaxis("vd"))
all_results_auto = all_results_manual.C
all_error_auto = all_results_manual.C
all_error_manual = all_results_manual.C
with figlist_var() as fl:
    for j in range(n_repeats):
        data = clean_data.C
        data.add_noise(fake_data_noise_std)
        # at this point, the fake data has been generated
        data.ft(["ph1", "ph2"], unitary=True)
        dt = diff(data.getaxis("t2")[r_[0, 1]]).item()
        data.ft("t2", shift=True)
        data /= sqrt(ndshape(data)["t2"]) * dt
        # note that frq_slice is re-determined for each repeat.  This is on
        # purpose.
        s_int, frq_slice = frequency_domain_integral(
            data,
            signal_pathway=signal_pathway,
            excluded_pathways=excluded_pathways,
            indirect="vd",
            fl=(fl if j == 0 else None),
            return_frq_slice=True,
        )
        all_results_auto["repeats", j] = s_int.data
        all_error_auto["repeats", j] = s_int.get_error()
        # {{{ manually calculate the error using manual bounds
        manual_integration = select_pathway(data, signal_pathway)[
            "t2":manual_bounds
        ].real
        if j == 0:  # the following are always the same for manual integration
            N = ndshape(manual_integration)["t2"]
            df = diff(data.getaxis("t2")[r_[0, 1]]).item()
        manual_integration.integrate("t2")
        all_results_manual["repeats", j] = manual_integration
        # Here, I use var, but I could also do this manually (w/out the ddof)
        # using abs( )**2 and then mean
        inactive_frq_datapoint_var = (
            select_pathway(data, {"ph1": 0, "ph2": 0})["t2":manual_bounds]
            .C.run(
                lambda x, axis=None: var(x, ddof=1, axis=axis) / 2, "t2"
            )  # the 2 is b/c var gives sum of real var and imag var
            .mean_all_but(["vd"])
        )
        all_error_manual["repeats", j] = (
            N * df**2 * inactive_frq_datapoint_var
        ).run(sqrt)
        # }}}
        print("#%d" % j)
    fl.next("different types of error")
    for this_label, this_data, this_error in [
        ("manual bounds", all_results_manual, all_error_manual),
        ("auto bounds", all_results_auto, all_error_auto),
    ]:
        this_data.mean("repeats", std=True)
        fl.plot(
            this_data,
            ".",
            capsize=6,
            label=this_label + " with true error",
            alpha=0.5,
        )
        # TODO ☐: for auto code, the propagated error is smaller than the true error
        this_error.mean("repeats")
        this_data.set_error(this_error.data)
        fl.plot(
            this_data,
            ".",
            capsize=6,
            label=this_label + " with propagated error",
            alpha=0.5,
        )
