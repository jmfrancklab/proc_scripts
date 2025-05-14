"""
Weighted Integral Function Example
==================================

The weighted integral function is utilized in integrating a signal in the time
domain.  Specifically here, we generate a sinc function whose FT is a heaviside
hat function with a width equal to the integration bounds and integrate in the
time domain.  These integrals are then compared to the same integrals
calculated in the frequency domain.
"""

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
    axis_coords=OrderedDict(
        [
            ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
            ("t2", nddata(r_[-npts:npts] * 30 / npts, "t2")),
        ]
    ),
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
    # {{{ Zero fill 2 x and make real/symmetric
    nPoints = len(select_pathway(data, signal_pathway)["repeats", 0].data)
    data.ft("t2", pad=2 * nPoints)
    data.run(np.real)
    data.ft_new_startpoint("t2", "t").set_ft_prop("t2", None)
    data.ift("t2", shift=True)
    # }}}
    data.ft("t2")
    # {{{ Normalize to integral of 1
    data /= (
        select_pathway(data["t2":frq_slice], signal_pathway)
        .C.real.integrate("t2")
        .mean("repeats")
        .item()
    )
    # }}}
    data.reorder(["ph1", "repeats", "t2"])
    for j in range(len(data["repeats"])):
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
    # {{{ Normalizing the sinc for integration
    for_norm = (mysinc.C**2).integrate("t2").item()
    mysinc /= np.sqrt(for_norm)
    mysinc *= np.sqrt(int_width)
    # }}}
    data.ift("t2")
    mysinc.set_units("t2", data.get_units("t2"))
    data.ft("t2").set_plot_color("tab:blue")
    mysinc.ft("t2").set_plot_color("tab:orange")
    fl.next("Frequency domain")
    fl.plot(select_pathway(data["repeats", 0], signal_pathway), alpha=0.5)
    fl.plot(mysinc * abs(data["repeats", 0]).max().item(), alpha=0.5)
    # }}}
    data.ift("t2")
    mysinc.ift("t2")
    fl.next("Time domain")
    fl.plot(select_pathway(data["repeats", 0], signal_pathway), alpha=0.5)
    fl.plot(mysinc / int_width, alpha=0.5)
    # {{{ Integrate in time domain
    apo_data = data * mysinc
    for j in range(len(data["repeats"])):
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
