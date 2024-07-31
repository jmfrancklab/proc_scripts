r"""
Calculate actual beta as a function of pulse length
===================================================
Assuming the data is the capture of the pulse sequence as seen on the GDS
oscilloscope (acquired using FLInst/examples/calib_pulses.py), 
here the data is converted to analytic power, frequency filtered
and the absolute is taken prior to integrating to return the beta where
:math:`\beta = \int \sqrt{P(t)} dt` 
"""
import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

color_cycle = cycle(
    [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
)

atten_ratio = 101.35  # attenutation ratio
show_all = False  # diagnostic
with psd.figlist_var() as fl:
    for filename, nodename, label in [
        ("240730_test_amp1_fin_pulse_calib.h5", "pulse_calib_1", "amplitude = 1"),
        ("240730_test_amp0p1_fin_pulse_calib.h5", "pulse_calib_1", "amplitude = 0.1"),
    ]:
        d = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        thiscolor = next(color_cycle)
        # {{{ allocate shapes for the final t_v_beta and beta_v_t plots
        t_v_beta = d.shape.pop("t").alloc().rename("p_90", "desired_beta")
        t_v_beta.setaxis("desired_beta", d.get_prop("desired_betas"))
        beta_v_t = d.shape.pop("t").alloc().rename("p_90", "At_p")
        beta_v_t.setaxis(
            "At_p",
            list(
                d.get_prop("acq_params")["amplitude"] * d.get_prop("set_p90s")[j]
                for j in range(len(d["p_90"]))
            ),
        )
        # }}}
        for j in range(len(d["p_90"])):
            s = d["p_90", j].C
            s.set_units("t", "s")
            s *= atten_ratio
            s /= np.sqrt(2)  # Vrms
            s /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)
            if show_all:
                fl.basename = r"pulse length = %s" % str(s.get_prop("set_p90s")[j])
                fl.next(fl.basename)
                # {{{ Plot raw capture
                fl.plot(s, alpha=0.2, color="blue")
                # }}}
            # {{{ make analytic
            s.ft("t", shift=True)
            s = s["t":(0, None)]
            s *= 2
            s["t":0] *= 0.5
            # }}}
            if show_all:
                s.ift("t")
                # {{{ Plot analytic
                fl.plot(abs(s), color="orange")
                s.ft("t")
                # }}}
            # {{{ frequency filter
            s["t":(0, 11e6)] *= 0
            s["t":(24e6, None)] *= 0
            # }}}
            s.ift("t")
            if show_all:
                # {{{ plot and integrate 90 pulse
                fl.plot(abs(s), color="red")
                # }}}
            int_range = abs(s).contiguous(lambda x: x > 0.01 * s.max())[0]
            # slightly expand int range to include rising edges
            int_range[0] -= 0.5e-6
            int_range[-1] += 0.5e-6
            beta1 = abs(s["t":int_range]).integrate("t").data.item() * 1e6
            beta_v_t["At_p", j] = beta1
            if show_all:
                plt.ylabel(r"$\sqrt{P_{pulse}}$")
                plt.text(
                    int_range[0] * 1e6 - 1,
                    -1,
                    r"$t_{90} \sqrt{P_{tx}} = %f s \sqrt{W}$" % beta1,
                )
                plt.axvline(int_range[0] * 1e6, ls=":", alpha=0.2)
                plt.axvline(int_range[1] * 1e6, ls=":", alpha=0.2)
        # {{{ make nddata for beta of 90 pulse and 180 pulse
        # {{{ beta vs t
        fl.next(r"Measured $\beta$ vs A * $t_{pulse}$")
        fl.plot(beta_v_t, "o", color=thiscolor, label=label)
        if beta_v_t < 6:
            mask = np.ones_like(beta_v_t.data, dtype=bool)
        else:
            mask = beta_v_t > 6.5
        calibration_data = psd.nddata(beta_v_t.data[mask], [-1], ["At_p"]).setaxis(
            "At_p", beta_v_t.getaxis("At_p")[mask]
        )
        calibration_data.sort("At_p")
        if beta_v_t < 6:
            c = calibration_data.polyfit("At_p", order=10)
        else:
            c = calibration_data.polyfit("At_p", order=1)
        fit_beta_v_t = np.polyval(c[::-1], beta_v_t.getaxis("At_p"))
        fit = psd.nddata(fit_beta_v_t, "At_p").setaxis("At_p", beta_v_t.getaxis("At_p"))
        fl.plot(fit, color=thiscolor, ls=":", alpha=0.5)
        psd.gridandtick(plt.gca())
        plt.ylabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.xlabel(r"programmed amplitude*$t_{pulse}$ / $\mu$s")
        # }}}
        # {{{ t vs beta
        t_v_beta = beta_v_t.shape.alloc().rename("At_p", "beta")
        t_v_beta.setaxis("beta", beta_v_t.data)
        t_v_beta.data[:] = beta_v_t.getaxis("At_p")
        fl.next(r"Amplitude*$t_{pulse}$ vs Measured $\beta$")
        fl.plot(t_v_beta, "o", color=thiscolor, label=label)
        if t_v_beta < 40:
            mask = np.ones_like(t_v_beta.data, dtype=bool)
        else:
            mask = t_v_beta > 45
        calibration_data = psd.nddata(t_v_beta.data[mask], [-1], ["beta"]).setaxis(
            "beta", t_v_beta.getaxis("beta")[mask]
        )
        calibration_data.sort("beta")
        if t_v_beta < 40:
            c = calibration_data.polyfit("beta", order=10)
        else:
            c = calibration_data.polyfit("beta", order=1)
        fit_t_v_beta = np.polyval(c[::-1], t_v_beta.getaxis("beta"))
        fit = psd.nddata(fit_t_v_beta, "beta").setaxis("beta", t_v_beta.getaxis("beta"))
        fl.plot(fit, color=thiscolor, ls=":", alpha=0.5)
        plt.xlabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.ylabel(r"programmed amplitude*$t_{pulse}$ / $\mu$s")
        psd.gridandtick(plt.gca())
        # }}}
