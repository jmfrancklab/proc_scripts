"""
Weighted Integral Function Example
==================================

The weighted integral function is utilized in integrating a signal.
Specifically here, we generate a sinc function whose FT is a heaviside hat
function with a width equal to the integration bounds.
"""
# TODO ☐: pulling very heavily from the stuff in your paper, can we
#         please just show a repeated series of simulated data, and
#         integrate both in the frequency domain and in the time domain.
#         Do this without apo, and use the same (manual) frequency range for both,
#         so that the results match.

from pyspecdata import r_, nddata, figlist_var, ndshape,fake_data
from pyspecProcScripts import select_pathway,apod_matched_filter, heaviside_time_domain, frequency_domain_integral
import numpy as np
from pylab import rcParams
from numpy.random import seed
import matplotlib.pyplot as plt
import sympy as s
from collections import OrderedDict

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1
excluded_pathways = [{"ph1":0},{"ph1":1},{"ph1":3}]
t2, repeats, ph1 = s.symbols("t2 repeats ph1")
signal_pathway = {"ph1":1}
n_repeats = 50
f_range = (-20,20)
int_width = f_range[1]-f_range[0]
# {{{ Generate data
data = fake_data(
        expression=(
            100*s.exp(+1j*2*s.pi*t2-t2*16*s.pi)
            +1e-10*repeats
            ),
        axis_coords=OrderedDict(
            [
                ("ph1",nddata(r_[0:4] / 4.0,"ph1")),
                ("t2", nddata(r_[0:1.024:1024j],"t2")),
                ("repeats",nddata(r_[0:n_repeats] + 1.0,"repeats")),
                ]
            ),
        signal_pathway = signal_pathway,
        fake_data_noise_std = 7.0,
        scale=10.0
        )
# }}}
# {{{ Note that we start at zero, but still need the echolike → causal
#     conversion
data = data["t2":(0,None)]
data *= 2
data["t2":0] *= 0.5
# }}}
t_stop_orig = data.getaxis("t2")[-1] # this will be used later in
#                                      making the cleaned up non zero
#                                      data
orig_causal_data = data.C # needed later for integration
# {{{ Zero fill 2 x and make real/symmetric
data.ft("t2",pad = 2*1024)
data.run(np.real)
data.ft_new_startpoint("t2","t").set_ft_prop("t2",None)
data.ift("t2",shift=True)
# }}}
data.ft("t2")
# {{{ Normalize to integral of 1
data /= (
        select_pathway(data["t2":f_range],signal_pathway)
        .real.integrate("t2")
        .mean("repeats")
        .item()
        )
# }}}
mysinc = heaviside_time_domain(select_pathway(data,signal_pathway),f_range)
# {{{ Normalizing the sinc for integration
mysinc.ft("t2")
for_norm = (mysinc.C**2).integrate("t2").item()
mysinc /= np.sqrt(for_norm)
mysinc.ift("t2")
mysinc *= np.sqrt(int_width)
# }}}
data.ift("t2")
mysinc.set_units("t2",data.get_units("t2"))
nonzero_data = data["t2":(-t_stop_orig,t_stop_orig)]
with figlist_var() as fl:
    data.ft("t2")
    data.set_plot_color("tab:blue")
    mysinc.ft("t2")
    mysinc.set_plot_color("tab:orange")
    fl.next("Frequency domain")
    fl.plot(select_pathway(data["repeats",0],signal_pathway), alpha = 0.5)
    fl.plot(mysinc*abs(data).mean("repeats").max().item(),alpha = 0.5)
    data.ift("t2")
    mysinc.ift("t2")
    fl.next("Time domain")
    fl.plot(select_pathway(data["repeats",0],signal_pathway),alpha = 0.5)
    fl.plot(mysinc/int_width,alpha = 0.5)
    # {{{ Integrate in time domain
    apo_data = data * mysinc
    # {{{ allocate array to place integrals in
    t_integral = ndshape(
        orig_causal_data.fromaxis("repeats")
    ).alloc()
    t_integral.setaxis("repeats",orig_causal_data["repeats"])
    # }}}
    for j in range(len(data["repeats"])):
        t_integral["repeats",j] = apo_data["repeats",j].C.real.integrate("t2").item()
    # }}}
    # {{{ Integrate in frequency domain
    data.ft("t2")
    f_integral, returned_frq_slice = frequency_domain_integral(
            data,
            signal_pathway = signal_pathway,
            excluded_frqs=[f_range],
            excluded_pathways = excluded_pathways,
            indirect="repeats",
            return_frq_slice = True)
    print(f_integral)
    print(t_integral)
    quit()

    # {{{ Plot
    for thisdata in [data, wt_fn]:
        fl.plot(thisdata, ax=thisax[0])
        thisdata.ift("t2")
        # Just normalize by the max so they are scaled correctly
        fl.plot(thisdata["t2":(-0.6, 0.6)] / thisdata.max(), ax=thisax[1])
    fl.show();quit()
    # }}}
