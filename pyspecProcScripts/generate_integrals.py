from pylab import *
from pyspecdata import *
import numpy as np
import matplotlib.pyplot as plt
from pyspecProcScripts import *
from .simple_functions import select_pathway
from .correlation_alignment import correl_align
from .integral_w_error import integral_w_errors
from .envelope import fit_envelope,L2G
from .DCCT_func import DCCT
import matplotlib.patches as patches

this_figsize = (6, 12)

def generate_integrals(
    s,
    searchstr="",
    signal_pathway={"ph1": 0, "ph2": 1},
    excluded_pathways=[(0, 0)],
    f_range=(None, None),
    direct="t2",
    indirect="indirect",
    slice_mult=2,
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
    slice_multiplier:   int
                        in determining the autoslice this is a multiplier where the higher
                        the value the wider the slice
    indirect:       str
                    Indirect dimension of the signal. Usually 'power' or 'nScans'
                    for example.
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
    raw_s = s
    # {{{DC offset correction
    # }}}
    s = fid_from_echo(s,signal_pathway,slice_multiplier = slice_mult, fl = fl)
    fl.next('FID f domain')
    fl.image(s)
    s.ift(direct)
    if clock_correction:
        # {{{clock correction
        clock_corr = nddata(np.linspace(-2, 2, 2500), "clock_corr")
        s.ft(direct)
        if indirect == 'vd':
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
    fl.next('FID t domain')
    fl.image(s)
    #get sigma for correlation alignment
    lambda_L = fit_envelope(s)
    s.ft('t2')
    # {{{ roughly align
    mysgn = determine_sign(select_pathway(s,signal_pathway))
    matched = (select_pathway(s,signal_pathway)*mysgn).ift('t2')
    matched *= exp(-pi*lambda_L*matched.fromaxis('t2'))
    matched.ft('t2')
    frq_atmax = matched.real.argmax('t2')
    s.ift('t2')
    t2 = s.fromaxis('t2')
    s *= exp(-1j*2*pi*frq_atmax*t2)
    s.ft('t2')
    fl.next('simple shift')
    fl.image(s)
    # }}}
    # {{{apply correlation alignment
    s.ift(list(signal_pathway))
    opt_shift,sigma,my_mask = correl_align(
            s*mysgn, indirect_dim=indirect,
            sigma = lambda_L,
            signal_pathway=signal_pathway)
    s.ift(direct)
    s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis(direct))
    s.ft(list(signal_pathway))
    s.ft(direct)
    fl.next('second alignment')
    fl.image(s)
    # }}}
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
                s.C.mean("nScans").real,
                signal_pathway,
                error_path,
                convolve_method="Gaussian",
                indirect=indirect,
                return_frq_slice=True,
            )
        else:
            s_int, frq_slice = integral_w_errors(
                s.real,
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
            x = s.getaxis(direct)
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
                select_pathway(s.real.run(complex128), signal_pathway),
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
    if indirect == "vd":
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
