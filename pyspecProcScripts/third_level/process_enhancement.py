"""
Processes enhancement data
==========================
Processes data acquired from an enhancement experiment 
and plots the resulting enhancement curve normalized.
"""
import pyspecdata as psd
from .phasing import hermitian_function_test
from .correlation_alignment import correl_align
from .calc_error import frequency_domain_integral
from sympy import symbols
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from .simple_functions import select_pathway

plt.rcParams.update(
    {
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
    }
)
t2 = symbols("t2")
thesecolors = cycle(list("bgrcmykw"))


def as_scan_nbr(s):
    return s.C.setaxis("power", "#").set_units("power", "scan #")


# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
def process_enhancement(
    s,
    searchstr="",
    signal_pathway={"ph1": 1},
    excluded_pathways=[(0, 0)],
    freq_range=(None, None),
    t_range=(0, 0.083),
    flip=False,
    sign=None,
    fl=None,
):
    s *= sign
    if fl is not None:
        fl.side_by_side(
            "show frequency limits\n$\\rightarrow$ use to adjust freq range",
            s,
            thisrange=freq_range,
        )  # visualize the frequency limits
    s.ift("t2")
    s.reorder(["ph1", "power", "t2"])
    if fl is not None:
        fl.push_marker()
        fl.next("time domain")
        fl.image(as_scan_nbr(s))
    plt.rcParams.update(
        {
            "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
            "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
            "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
        }
    )
    s.ift(["ph1"])
    # {{{ Applying DC offset correction
    t_start = t_range[-1] / 4
    t_start *= 3
    rx_offset_corr = s["t2":(t_start, None)]
    rx_offset_corr = rx_offset_corr.data.mean()
    s -= rx_offset_corr
    s.ft("t2")
    s.ft(["ph1"])
    # }}}
    zero_crossing = (
        abs(select_pathway(s, signal_pathway))
        .sum("t2")
        .argmin("power", raw_index=True)
        .item()
    )
    s = s["t2":freq_range]
    if fl is not None:
        fl.next("freq_domain before phasing")
        fl.image(s.C.setaxis("power", "#").set_units("power", "scan #"))
    # {{{Applying phasing corrections
    s.ift("t2")  # inverse fourier transform into time domain
    best_shift, max_shift = hermitian_function_test(
        select_pathway(s, signal_pathway).C.convolve("t2", 3e-4)
    )
    best_shift = 0.033e-3
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    psd.logger.info(psd.strm("applying zeroth order correction"))
    s.ift(["ph1"])
    phasing = s["t2", 0].C
    phasing.data *= 0
    phasing.ft(["ph1"])
    phasing["ph1", 1] = 1
    phasing.ift(["ph1"])
    s /= phasing
    ph0 = s["t2":0] / phasing
    ph0 /= abs(ph0)
    s /= ph0
    s.ft(["ph1"])
    psd.logger.info(psd.strm(s.dimlabels))
    s.ft("t2")
    if fl is not None:
        fl.next("phase corrected")
        fl.image(as_scan_nbr(s))
    s.reorder(["ph1", "power", "t2"])
    psd.logger.info(psd.strm("zero corssing at", zero_crossing))
    # }}}
    # {{{Correcting power axis
    # print(s.getaxis('power'))
    # quit()
    # power_axis_dBm = array(s.get_prop('meter_powers'))
    # print(power_axis_dBm)
    # power_axis_W = zeros_like(power_axis_dBm)
    # power_axis_W[:] = 10**(power_axis_dBm/10)
    # power_axis_W = r_[0,power_axis_W]
    # print(power_axis_W)
    # quit()
    # s.setaxis('power',power_axis_W)
    # s.set_units('power','W')
    # }}}
    # {{{Applying correlation alignment
    s.ift(["ph1"])
    opt_shift, sigma = correl_align(
        s, indirect_dim="power", ph1_selection=1, sigma=0.001
    )
    s.ift("t2")
    s *= np.exp(-1j * 2 * np.pi * opt_shift * s.fromaxis("t2"))
    s.ft("t2")
    fl.basename = None
    if fl is not None:
        fl.next(r"after correlation, $\varphi$ domain")
        fl.image(as_scan_nbr(s))
    s.ift("t2")
    s.ft(["ph1"])
    if fl is not None:
        fl.next("after correlation alignment FTed ph")
        fl.image(as_scan_nbr(s))
    s.reorder(["ph1", "power", "t2"])
    if fl is not None:
        fl.next("after correlation -- time domain")
        fl.image(as_scan_nbr(s))
    s.ft("t2")
    if fl is not None:
        fl.next("after correlation -- frequency domain")
        fl.image(as_scan_nbr(s))
    # }}}
    s.ift("t2")
    d = s.C
    d.ft("t2")
    d.ift("t2")
    d = d["t2" : (0, t_range[-1])]
    d["t2":0] *= 0.5
    d.ft("t2")
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a
    #     compact method for N pulses
    error_pathway = (
        set(((j) for j in range(psd.ndshape(d)["ph1"])))
        - set(excluded_pathways)
        - set([(signal_pathway["ph1"])])
    )
    error_pathway = [{"ph1": j} for j in error_pathway]
    # }}}
    # {{{ integrating with error bar calculation
    d_, frq_slice, std = frequency_domain_integral(
        d,
        signal_pathway,
        error_pathway,
        indirect="power",
        fl=fl,
        return_frq_slice=True,
    )
    x = d_.get_error()
    x[:] /= np.sqrt(2)
    d = d_.C
    # }}}
    # {{{Normalizing by max
    idx_maxpower = np.argmax(s.getaxis("power"))
    d /= max(d.data)
    # }}}
    power_axis_dBm = np.array(s.get_prop("meter_powers"))
    power_axis_W = np.zeros_like(power_axis_dBm)
    power_axis_W[:] = 1e-2 * 10 ** ((power_axis_dBm[:] + 10.0) * 1e-1)
    power_axis_W = psd.r_[0, power_axis_W]
    d.setaxis("power", power_axis_W)
    # d.set_units('power','W')
    if flip:
        d = 1 - d
    if fl is not None:
        fl.next("E(p)")
        fl.plot(d["power", :-3], "ko", capsize=6, alpha=0.3)
        fl.plot(d["power", -3:], "ro", capsize=6, alpha=0.3)
        fl.pop_marker()
    enhancement = d
    return enhancement, idx_maxpower
