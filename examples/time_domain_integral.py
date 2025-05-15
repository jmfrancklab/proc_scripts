"""
Time-domain integral function
=============================

The weighted integral function is utilized in integrating a signal in the time
domain.  Specifically here, we generate a sinc function whose FT is a heaviside
hat function with a width equal to the integration bounds and integrate in the
time domain.  These integrals are then compared to the same integrals
calculated in the frequency domain.
"""

# TODO ☐: the wording above is weird to me -- I don't understand why you make
#         this about a "weighted integral function".  The point is to show that
#         your time-domain integration matches the frequency-domain
#         integration.

from pyspecdata import r_, nddata, figlist_var, ndshape, fake_data
from pyspecProcScripts import (
    select_pathway,
    heaviside_time_domain,
)
import numpy as np
from pylab import rcParams
from numpy.random import seed
import matplotlib.pyplot as plt
import sympy as s
from collections import OrderedDict

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1
t2, repeats, ph1 = s.symbols("t2 repeats ph1")
signal_pathway = {"ph1": 1}
n_repeats = 100
npts = 8192
frq_slice = (-30, 30)
int_width = frq_slice[1] - frq_slice[0]
# {{{ Generate data with an amplitude of 21 and t constant = 320 ms
d = fake_data(
    expression=(21 * s.exp(+1j * 2 * s.pi * t2 - abs(t2) * 0.32 * s.pi)),
    axis_coords=OrderedDict([
        ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
        ("t2", nddata(r_[-npts:npts] * 30 / npts, "t2")),
    ]),
    signal_pathway=signal_pathway,
    scale=0,
    fake_data_noise_std=0,
)
data = (
    ndshape([len(d["t2"]), len(d["ph1"]), n_repeats], ["t2", "ph1", "repeats"])
    .alloc()
    .setaxis("t2", d.getaxis("t2"))
    .setaxis("ph1", d.getaxis("ph1"))
    .set_units("t2", "s")
    .setaxis("repeats", r_[0:n_repeats])
)
for j in range(len(data.getaxis("repeats"))):
    data["repeats", j] = d.C
    data["repeats", j].add_noise(0.3)
# }}}
# {{{ Allocate an nddata for the t_integrals
t_results = (
    ndshape([n_repeats], ["repeats"])
    .alloc()
    .setaxis("repeats", r_[1 : n_repeats + 1])
)
f_results = t_results.C
# }}}
# {{{ Center data at 0 Hz
data.ft("t2", shift=True)
nu_0 = (
    select_pathway(data, signal_pathway).C.mean("repeats").argmax("t2").item()
)
data.reorder(["repeats", "t2"])
data.ift("t2")
data *= np.exp(-1j * 2 * np.pi * nu_0 * data.fromaxis("t2"))
# }}}
# {{{ Note that we start at zero, but still need the echolike → causal
#     conversion
data = data["t2":(0, None)]
data *= 2
data["t2":0] *= 0.5
# }}}
with figlist_var() as fl:
    # TODO ☐: The goal is to (1) take an integral in the frequency domain (2)
    #         take an integral in the time domain (3) show that the integrals
    #         and the errors match.  This is currently showing that the
    #         integrals match, but there are no errors.  Given that you are
    #         using a fixed bounds (good, keeps it simple), I would just call
    #         the function that calculates the masked variance, and calculate
    #         the error from that.
    # {{{ Zero fill 2 x and make real/symmetric
    nPoints = len(select_pathway(data, signal_pathway)["repeats", 0].data)
    data.ft("t2", pad=2 * nPoints)
    data.run(np.real)
    data.ft_new_startpoint("t2", "t").set_ft_prop("t2", None).ift(
        "t2", shift=True
    )
    # }}}
    # TODO ☐: I removed the normalization b/c you choose the size of the data
    #         above.  Get it to work w/out rolling this back.
    data.ft("t2")
    data.reorder(["ph1", "repeats", "t2"])
    for j in range(len(data["repeats"])):
        # TODO ☐: why is there a loop here? Vectorize your code!
        f_results["repeats", j] = (
            select_pathway(data["repeats", j]["t2":frq_slice], signal_pathway)
            .C.real.integrate("t2")
            .item()
        )
    # {{{ Make sinc function that will be used as the weighted integral
    #     function
    mysinc = heaviside_time_domain(
        select_pathway(data, signal_pathway), frq_slice
    )
    # TODO ☐: the extra plots were unneeded and we don't want to add a weird
    #         extra normalization to our w(t), so I deleted those.  Get code to
    #         work without rolling back that change.
    # }}}
    data.ift("t2")
    # TODO ☐: the following should not be needed here -- you might want to edit heaviside_time_domain
    mysinc.set_units("t2", data.get_units("t2"))
    # {{{ Integrate in time domain
    # TODO ☐: objecting to your variable naming here -- this is NOT apodized!
    #         You are just calculating the integral.  Aside from that, why not
    #         just do the expression in the loop (which should not be a loop!)
    #         together with the multiplication in one expression.  That matches
    #         your math better, anyways.
    apo_data = data * mysinc
    for j in range(len(data["repeats"])):
        # TODO ☐: why is there a loop here? Vectorize your code!
        t_results["repeats", j] = (
            select_pathway(apo_data["repeats", j], signal_pathway)
            .C.real.integrate("t2")
            .item()
        )
    # }}}
    avg_t_integral = np.mean(t_results.data).real.item()
    avg_f_integral = np.mean(f_results.data).real.item()
    print(
        "The integral of the first repeat in the frequency domain is %.5f"
        " which is approximately equal to the integral of the first repeat in"
        " the time domain which is %.5f"
        % (
            f_results["repeats", 0].real.item(),
            t_results["repeats", 0].real.item(),
        )
    )
    print(
        "Meanwhile the average integral of all repeats taken in the frequency"
        " domain is %.5f and when taken in the time domain the average is %.5f"
        % (
            avg_t_integral,
            avg_f_integral,
        )
    )
    # {{{ Plot integrals from frequency domain integration and from time domain
    #     integration
    fl.next("Integrals")
    fl.plot(t_results, "o", color="tab:blue", label="time domain integration")
    fl.plot(
        f_results.set_error(None),
        "o",
        color="tab:orange",
        label="frequency domain integration",
    )
    plt.axhline(avg_t_integral, color="tab:blue")
    plt.axhline(avg_f_integral, color="tab:orange")
    # }}}
