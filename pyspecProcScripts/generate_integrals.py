from pylab import *
from pyspecdata import *
from scipy.optimize import minimize, leastsq
from sympy import exp as s_exp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, latex, Symbol
from pyspecProcScripts import *
from .simple_functions import select_pathway
from .correlation_alignment import correl_align
from .integral_w_error import integral_w_errors
from .DCCT_func import DCCT
import matplotlib.patches as patches

this_figsize = (6, 12)
t2 = symbols("t2")


def generate_integrals(
    s,
    searchstr="",
    signal_pathway={"ph1": 0, "ph2": 1},
    excluded_pathways=[(0, 0)],
    f_range=(None, None),
    direct="t2",
    indirect="indirect",
    alias_slop=3,
    clock_correction=True,
    error_bars=True,
    correlate=True,
    fl=None,
):
    """Applies appropriate phase corrections, along with DC offset corrections,
    alignment, and integrates the data to transform 2D datasets into a table
    of integrals.

    Parameters
    ==========
    searchstr:      str
                    string for the title of the plots produced.
    signal_pathway: dict
                    Dictionary of the path of the desired signal.
    excluded_pathways:  dict
                        Dictionary of all coherence pathways that are not
                        the signal pathway.
    f_range:        tuple
                    Frequency range over which the signal resides.
    direct:         str
                    Direct dimension of the signal
    indirect:       str
                    Indirect dimension of the signal. Usually 'power' or 'nScans'
                    for example.
    alias_slop:     int
                    Aliasing_slop used in the hermitian function.
    clock_correction:   bool
                        If true, will apply a clock correction. Is needed for IR
                        but not enhancement data.
    error_bars:     bool
                    Option to apply integration *with* error bars.
    correlate:      bool
                    Option to apply correlation alignment to data.

    Returns
    =======
    s:              nddata
                    Data with applied corrections but *not* integrated
    s_int:          nddata
                    Integrated and corrected data
    """
    s.ift(direct)
    s.reorder([indirect, direct], first=False)
    # {{{DC offset correction
    s.ift(
        list(signal_pathway)
    )  # list(signal_pathway) takes the keys of the signal_pathway dictionary
    t_max = s.getaxis("t2").max()
    rx_offset_corr = s["t2" : (t_max * 0.75, None)]
    rx_offset_corr = rx_offset_corr.mean(["t2"])
    s -= rx_offset_corr
    s.ft(direct)
    s.ft(list(signal_pathway))
    # }}}
    # {{{phase correction
    s = s[direct:f_range]
    s.ift(list(signal_pathway))
    raw_s = s.C  # will be used for imaging raw data with proper scaling
    s.ft(list(signal_pathway))
    s.ift(direct)
    best_shift = hermitian_function_test(
        select_pathway(s, signal_pathway), aliasing_slop=alias_slop
    )
    logger.info(strm("best shift is", best_shift))
    s.setaxis(direct, lambda x: x - best_shift).register_axis({direct: 0})
    s.ft(direct)
    if "nScans" in s.dimlabels:
        s_mean = s.C.mean("nScans")
    else:
        s_mean = s
    mysgn = determine_sign(select_pathway(s_mean, signal_pathway))  # must be
    # determined after hermitian phasing, or else you can actually get a
    # sign flip when the signal drifts far enough in frequency space
    s /= zeroth_order_ph(select_pathway(s, signal_pathway))
    ph_corr_s = (
        s.C
    )  # will be used for imaging phased data with proper scaling later
    s.ift(direct)
    # }}}
    # {{{Correlate
    if correlate:
        s.ft(direct)
        s.ift(list(signal_pathway))
        opt_shift, sigma, my_mask = correl_align(
            s * mysgn, indirect_dim=indirect, signal_pathway=signal_pathway
        )
        s.ift(direct)
        s *= np.exp(-1j * 2 * pi * opt_shift * s.fromaxis(direct))
        s.ft(list(signal_pathway))
        s.ft(direct)
        scale_factor = abs(s.C).max().item()
    # }}}
    if fl:
        this_fig = figure(figsize=(20, 10))
        DCCT(
            raw_s,
            this_fig,
            total_spacing=0.1,
            RHS_pad=0.76,
            allow_for_text_default=5,
            allow_for_ticks_default=50,
            text_height=50,
            custom_scaling=True,
            scaling_factor=scale_factor,
            plot_title="Raw Data \n for %s" % searchstr,
        )
        DCCT(
            ph_corr_s,
            this_fig,
            total_spacing=0.1,
            LHS_pad=0.25,
            RHS_pad=0.51,
            allow_for_text_default=5,
            allow_for_ticks_default=50,
            text_height=50,
            custom_scaling=True,
            scaling_factor=scale_factor,
            plot_title="Phased Data \nfor %s" % searchstr,
        )
        if correlate:
            DCCT(
                s,
                this_fig,
                total_spacing=0.1,
                custom_scaling=True,
                allow_for_text_default=5,
                allow_for_ticks_default=50,
                LHS_pad=0.5,
                RHS_pad=0.27,
                scaling_factor=scale_factor,
                text_height=50,
                plot_title="Aligned Data \nfor %s" % searchstr,
            )
    s.ift(direct)
    if clock_correction:
        # {{{clock correction
        clock_corr = nddata(np.linspace(-2, 2, 2500), "clock_corr")
        s.ft(direct)
        if "vd" is indirect:
            s_clock = s["ph1", 0]["ph2", 1].sum(direct)
        else:
            s_clock = s["ph1", 0].sum(direct)
        s.ift(list(signal_pathway))
        min_index = abs(s_clock).argmin(indirect, raw_index=True).item()
        s_clock *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
        s_clock[indirect, : min_index + 1] *= -1
        s_clock.sum(indirect).run(abs)
        clock_corr = s_clock.argmax("clock_corr").item()
        s *= np.exp(-1j * clock_corr * s.fromaxis(indirect))
        s.ft(list(signal_pathway))
        s.ift(direct)
        # }}}
    s_after = s.C
    s_after.ft(direct)
    s_after.ift(direct)
    s_after = s_after[direct:(0, None)]
    s_after[direct, 0] *= 0.5
    s_after.ft(direct)
    if fl:
        DCCT(
            s_after,
            this_fig,
            custom_scaling=True,
            total_spacing=0.2,
            allow_for_text_default=5,
            allow_for_ticks_default=50,
            LHS_pad=0.74,
            scaling_factor=scale_factor,
            text_height=50,
            plot_title="FID \nfor %s" % searchstr,
        )
    if "ph2" in s.dimlabels:
        logger.info(strm("PH2 IS PRESENT"))
        error_path = (
            set(
                (
                    (j, k)
                    for j in range(ndshape(s)["ph1"])
                    for k in range(ndshape(s)["ph2"])
                )
            )
            - set(excluded_pathways)
            - set([(signal_pathway["ph1"], signal_pathway["ph2"])])
        )
        error_path = [{"ph1": j, "ph2": k} for j, k in error_path]
    else:
        error_path = (
            set(((j) for j in range(ndshape(s)["ph1"])))
            - set(excluded_pathways)
            - set([(signal_pathway["ph1"])])
        )
        error_path = [{"ph1": j} for j in error_path]
        # {{{Integrate with error bars
    if error_bars:
        if "nScans" in s.dimlabels:
            s_int, frq_slice = integral_w_errors(
                s_after.C.mean("nScans"),
                signal_pathway,
                error_path,
                convolve_method="Gaussian",
                indirect=indirect,
                return_frq_slice=True,
            )
        else:
            s_int, frq_slice = integral_w_errors(
                s_after,
                signal_pathway,
                error_path,
                convolve_method="Gaussian",
                indirect=indirect,
                return_frq_slice=True,
            )
        x = s_int.get_error()
        x[:] /= sqrt(2)
        if fl is not None:
            fig0 = figure(figsize=this_figsize)
            x = s_after.getaxis(direct)
            dx = x[1] - x[0]
            fl.next("Real Integrated Data with Bounds", fig=fig0)
            ax1 = plt.axes([0.2933, 0.15, 0.6567, 0.713])
            yMajorLocator = lambda: mticker.MaxNLocator(
                nbins="auto", steps=[1, 2, 5, 10]
            )
            majorLocator = lambda: mticker.MaxNLocator(
                nbins="auto", steps=[1, 2, 2.5, 5, 10]
            )
            minorLocator = lambda: mticker.AutoMinorLocator(n=5)
            ax1.xaxis.set_major_locator(majorLocator())
            ax1.xaxis.set_minor_locator(minorLocator())
            ax1.set_ylabel(None)
            fl.image(
                select_pathway(s_after.real.run(complex128), signal_pathway),
                black=False,
                ax=ax1,
            )
            x1 = x[0]
            x2 = frq_slice[-1]
            wide = frq_slice[0] - x1 - dx
            wide2 = (x[-1] + dx) - x2
            y_bottom, y_top = ax1.get_ylim()
            p = patches.Rectangle(
                (x1, y_bottom),
                wide,
                y_top - y_bottom,
                angle=0.0,
                linewidth=1,
                fill=None,
                hatch="//",
                ec="k",
            )
            ax1.add_patch(p)
            q = patches.Rectangle(
                (x2, y_bottom),
                wide2,
                y_top - y_bottom,
                angle=0.0,
                linewidth=1,
                fill=None,
                hatch="//",
                ec="k",
            )
            ax1.add_patch(q)
    if indirect is "vd":
        s_int = s_int
        if fl is not None:
            fig1 = figure()
            fl.next("Integrated Data", fig=fig1)
            fl.plot(s_int, "o", capsize=6, alpha=0.3)
    else:
        s_int[indirect, :] /= s_int.data[0]
        if fl is not None:
            fig1 = figure()
            fl.next("Integrated Data", fig=fig1)
            fl.plot(s_int["power", :-3], "o", capsize=6, alpha=0.3)
    # }}}
    return s_int, s
