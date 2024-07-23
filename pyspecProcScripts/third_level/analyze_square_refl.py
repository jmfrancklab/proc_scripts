""" Analyze square wave reflection data
=======================================
Analyzes data acquired by applying a square wave pulse to an empty 
or filled NMR probe and acquiring the reflection profile.
"""

import pylab as plb
from pylab import r_, pi
import pyspecdata as psp
from scipy.signal import tukey
from scipy.optimize import minimize, leastsq
import sympy as s
import matplotlib.pyplot as plt
from pyspecdata import *


def analyze_square_refl(
    d,
    label="",
    fl=None,
    frq_bw=15e6,
    keep_after_pulse=2e-6,
    blip_range=[-0.1e-6, None],
    show_analytic_signal_phase=True,
    show_analytic_signal_real=False,
):
    r"""
    Plots the reflected square wave with a fit that calculates
    the Q for the probe being tested.

    Parameters
    ==========
    fl: figlist_var child class
        In addition to standard figlist_var methods, this must also include a
        `complex_plot` method that returns list "colors"
    d:      nddata
    label:  str
            appropriate name for dataset.
    frq_bw: int
            bandwidth in Hz.
    keep_after_pulse:   int
                        Amount of the time axis after the pulse
                        that you want to remain-changes depending
                        on the Q.
    blip_range:         lst
                        Range in which the signal blip exists.
    show_analytic_signal_phase: bool
                                Option on final plot to show the
                                phase of the signal.
    show_analytic_signal_real:  bool
                                Option on final plot to show the
                                real of the signal for comparison.
    """
    if len(label) > 0:
        fl.basename = label
    orig_basename = fl.basename
    if blip_range[-1] is None:
        blip_range[-1] = keep_after_pulse
    d.ft("t", shift=True)
    d.ift("t")
    d.ft("t")
    d = d["t":(-50e6, 50e6)]  # slice out a reasonable range
    d["t":(None, 0)]["t", :-1] = 0
    d *= 2  # multiply data by 2 because the equation
    #                       1/2a*exp(iwt)+aexp(-iwt) and the 2 negated the
    #                       half.
    if fl is not None:
        fl.next("show the frequency distribution")
        for j in d.getaxis("ch"):
            fl.plot(
                abs(d)["ch":j][lambda x: abs(x) > 1e-10],
                alpha=0.5,
                plottype="semilogy",
                label=f"CH{j} {label}",
            )
    frq_guess = abs(d["ch", 0]).argmax("t").item()  # find the peak
    frq_range = r_[-frq_bw, frq_bw] * 0.5 + frq_guess
    d["t" : (0, frq_range[0])] = 0  # set everything else to zero
    d["t" : (frq_range[1], None)] = 0
    # {{{ shouldn't have to do it this way, but something weird going on w/ aligndata
    tukey_filter = d.fromaxis("t").C
    tukey_filter = tukey_filter["t" : tuple(frq_range)].run(
        lambda x: tukey(len(x))
    )
    d["t" : tuple(frq_range)] *= tukey_filter
    if fl is not None:
        for j in d.getaxis("ch"):
            fl.plot(
                abs(d)["ch":j][lambda x: abs(x) > 1e-10],
                alpha=0.5,
                plottype="semilogy",
                label=f"CH{j} filtered {label}",
            )
        fl.grid()
    df = plb.diff(d.getaxis("t")[r_[0, 1]]).item()
    d.ift("t")
    # {{{ slice out pulse and surrounding and start time axis with pulse
    dt = plb.diff(d.getaxis("t")[r_[0, 1]]).item()
    pulse_slice = d["ch", 0].contiguous(
        lambda x: abs(x) > 0.5 * abs(x).data.max()
    )[
        0
    ]  # defines pulse slice based on control signal
    d.setaxis("t", lambda x: x - pulse_slice[0]).register_axis(
        {"t": 0}
    )  # set t=0 to the beginning of the pulse
    pulse_slice -= pulse_slice[0]
    d = d[
        "t" : tuple(pulse_slice + r_[-0.5e-6, keep_after_pulse])
    ]  # take a slice that starts 0.5 μs before the pulse and ends 2 μs after the pulse
    # }}}
    # {{{ determine ∂φ/∂t during the pulse, use to det frq
    pulse_middle_slice = tuple(pulse_slice + r_[+0.5e-6, -0.5e-6])
    pulse_middle = d["ch", 0]["t":pulse_middle_slice]
    ph_diff = pulse_middle["t", 1:] / pulse_middle["t", :-1]
    ph_diff *= abs(pulse_middle["t", :-1])
    # at this point, we have something of the form
    # ρe^{iΔφ} for each point -- average them
    ph_diff.sum("t")
    ph_diff = ph_diff.angle.item()  # Δφ above
    frq = ph_diff / dt / 2 / pi
    # }}}
    logger.info(strm("frq:", frq))
    d *= plb.exp(-1j * 2 * pi * frq * d.fromaxis("t"))  # mix down
    d.ft("t")
    if fl is not None:
        fl.next("after slice and mix down freq domain")
    if fl is not None:
        fl.plot(abs(d), label=label)
    d.ift("t")
    # {{{ zeroth order phase correction
    ph0 = d["ch", 0]["t":pulse_middle_slice].C.mean("t").item()
    pulse_middle_amp = abs(ph0)
    ph0 /= pulse_middle_amp
    d /= ph0
    # }}}
    if fl is not None:
        fl.basename = None
        fl.next("analytic signal")
        colors = fl.complex_plot(
            d,
            label=label,
            show_phase=show_analytic_signal_phase,
            show_real=show_analytic_signal_real,
            alpha=0.5,
        )
        # {{{ print the carrier
        # transform goes to "display", which is pixels
        # "inverted" goes back
        # here, I could set a mixed transformation,
        # but I think it's more understandable to just do manually
        ax = plb.gca()
        logger.info(strm("the amplitude is", pulse_middle_amp))
        _, y = ax.transData.transform(
            r_[0.0, pulse_middle_amp]
        )  # _ and y are now the display coords for the data coords of 0
        # and the middle of the pulse amp. same "spot" but now in the display
        # coordinate system.
        fontsize = 16
        nfigures = len(fl.figurelist)
        y -= fontsize * (nfigures / 3)  # adjusting y in display coordinates
        _, y = ax.transAxes.inverted().transform(
            r_[0, y]
        )  # going now from display coordinate system into Axes coordinate system
        plb.text(
            x=0.5,
            y=y,
            s=r"$\nu_{Tx}=%0.6f$ MHz" % (frq / 1e6),
            va="top",
            ha="center",
            size=fontsize,
            transform=ax.transAxes,  # display
            color=colors[0],
        )
        # }}}
        fl.basename = orig_basename
    scalar_refl = (
        d["ch", 1]["t" : (keep_after_pulse, pulse_middle_slice[-1])]
        .mean("t")
        .item()
    )
    if label == "hairpin probe":
        color = "darkorange"
    if label == "solenoid probe":
        color = "red"
    if fl is not None:
        fl.basename = None
        fl.next("blips")
        first_blip = (
            -d["ch", 1:2]["t" : tuple(blip_range)] + scalar_refl
        )  # correcting first blip
        fl.complex_plot(
            first_blip,
            "first" + "_" + label,
            show_phase=False,
            show_real=False,
            alpha=0.2,
            linestyle="--",
            linewidth=1,
            color=color,
        )
        secon_blip = d["ch", 1:2][
            "t" : tuple(blip_range + pulse_slice[1])
        ].setaxis("t", lambda x: x - pulse_slice[1])
        colors = fl.complex_plot(
            secon_blip,
            "second" + "_" + label,
            show_phase=False,
            show_real=False,
            alpha=0.8,
            color=color,
        )
    secon_blip = secon_blip[
        "ch", 0
    ]  # we need the ch axis for the complex plot,
    #                                  but it complicates things now
    decay = abs(secon_blip)
    decay_start = decay.C.argmax("t").item()
    decay = decay["t":(decay_start, None)]
    f = psp.fitdata(decay)
    A, R, C, t = s.symbols("A R C t", real=True)
    f.functional_form = A * s.exp(-t * R) + C
    f.set_guess({A: 0.3, R: 1 / 20 * 2 * pi * frq})
    f.fit()
    logger.info(strm("output:", f.output()))
    logger.info(strm("latex:", f.latex()))
    # {{{ calculate frequency offset
    decay_timescale = 3.0 / f.output("R")
    end_blip = (54.0 / f.output("R") * 2 * pi * frq) * 1e-9
    phases = secon_blip["t":(150e-9, end_blip)]
    # should actually determine error from noise std
    # -- the following is just an estimate
    phases.set_error(0.01)
    # should use new phdiff function in pySpecData to calculate the
    # following
    phase_diff = phases["t", 1:] / phases["t", :-1]
    dt = np.diff(phases.getaxis("t")[:2]).item()
    phase_diff.mean("t")  # favors points with a greater magnitude
    frq_offset = phase_diff.angle.item() / dt / 2 / pi
    logger.info(strm(frq_offset.real))
    frq_line_plot = phases.fromaxis("t")
    frq_line_plot -= frq_line_plot["t", 0]
    frq_line_plot *= 2 * pi * frq_offset  # so this contains 2πνt
    frq_line_plot = np.exp(1j * frq_line_plot)
    frq_line_plot *= phases["t", 0]
    if fl is not None:
        ax2 = fl.twinx(orig=False)
        fl.plot(
            phases.angle / 2 / pi,
            ".",
            linewidth=1,
            color=colors[-1],
            alpha=0.3,
            label="reflected angle " + label,
        )
        fl.plot(
            frq_line_plot.angle.set_error(None) / 2 / pi,
            color=colors[-1],
            linewidth=1,
        )
        ax2 = fl.twinx(orig=False)
        fl.plot(
            phases.angle / 2 / pi,
            ".",
            linewidth=1,
            color=colors[-1],
            alpha=0.3,
            label="reflected angle " + label,
        )
        fl.plot(
            frq_line_plot.angle.set_error(None) / 2 / pi,
            color=colors[-1],
            linewidth=1,
        )
        x = plb.mean(phases.getaxis("t")) / 1e-9
        y = np.angle(plb.mean(phases.data)) / 2 / pi
        ax2.text(
            x=x,
            y=y,
            s=" " * 5 + r"$\Delta\nu=%0.3g$ kHz" % (frq_offset / 1e3),
            va="bottom",
            ha="left",
            size=fontsize,
            transform=ax2.transData,
            color=color,
        )
        ax = fl.twinx(orig=True)
    # print('frq_offset',frq_offset,"for",label)
    # }}}
    Q = 1.0 / f.output("R") * 2 * pi * frq
    if fl is not None:
        ax = fl.twinx(orig=True)
        x = (54e-9 + 4.0 / f.output("R")) / 1e-9
        y = 0.1
        ax.text(
            x=x,
            y=y,
            s=" " * 5 + r"$Q=%d$" % (Q),
            va="bottom",
            ha="left",
            size=fontsize,
            transform=ax.transData,
            color=color,
        )
        fl.plot(
            f.eval(100).set_units("t", "s"), "k--", alpha=0.8
        )  # , label='fit, Q=%0.1f'%Q)
        fl.twinx(orig=True)
        fl.plot(f.eval(100).set_units("t", "s"), "k--", alpha=0.8)
        ax.grid(False)
        ax2.set_ylabel("phase / cyc", size=10)
        ax2.set_ylim(-0.5, 0.5)
        gridandtick(ax2)
