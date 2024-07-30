"""
Calculate actual beta as a function of pulse length
===================================================
Assuming the data is the capture of the pulse sequence as seen on the GDS
oscilloscope, here the data is converted to analytic power, frequency filtered
and the absolute is taken prior to integrating to return the beta value (us
sqrt(W)). 
"""
import pyspecdata as psd
import matplotlib.pyplot as plt
from numpy import sqrt

atten_ratio = 101.35  # attenutation ratio
show_all = False
with psd.figlist_var() as fl:
    for filename, nodename, label in [
        ("240730_test_amp1_fin_pulse_calib.h5", "pulse_calib_1", "amplitude = 1"),
        ("240730_test_amp0p1_fin_pulse_calib.h5", "pulse_calib_1", "amplitude = 0.1"),
    ]:
        d = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        beta = psd.ndshape(d).pop("t").alloc().rename("p_90", "prog_beta")
        beta.setaxis("prog_beta", d.get_prop("desired_betas"))
        t_v_beta = psd.ndshape(d).pop("t").alloc().rename("p_90", "At_p")
        t_v_beta.setaxis(
            "At_p",
            list(
                d.get_prop("acq_params")["amplitude"] * d.get_prop("set_p90s")[j]
                for j in range(len(d["p_90"]))
            ),
        )
        for j in range(len(d["p_90"])):
            s = d["p_90", j].C
            s.set_units("t", "s")
            s *= atten_ratio
            s /= sqrt(2)  # Vrms
            s /= sqrt(50)  # V/sqrt(R) = sqrt(P)
            if show_all:
                fl.basename = r"programmed $\beta$ = %s" % str(
                    s.get_prop("set_betas")[j]
                )
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
            beta90_int = abs(s).contiguous(lambda x: x > 0.01 * s.max())[0]
            beta1 = abs(s["t":beta90_int]).integrate("t").data.item() * 1e6
            beta["prog_beta", j] = beta1
            t_v_beta["At_p", j] = beta1
            if show_all:
                plt.ylabel(r"$\sqrt{P_{pulse}}$")
                plt.text(
                    beta90_int[0] * 1e6 - 1,
                    -1,
                    r"$t_{90} \sqrt{P_{tx}} = %f s \sqrt{W}$" % beta1,
                )
                plt.axvline(beta90_int[0] * 1e6, ls=":", alpha=0.2)
                plt.axvline(beta90_int[1] * 1e6, ls=":", alpha=0.2)
        # {{{ make nddata for beta of 90 pulse and 180 pulse
        fl.next(r"Measured $\beta$ vs programmed $\beta$")
        fl.plot(beta, "o", label=label)
        psd.gridandtick(plt.gca())
        plt.ylabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.xlabel(r"programmed $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        fl.next(r"Measured $\beta$ vs Amplitude*$t_{plse}$")
        fl.plot(t_v_beta, "o", label=label)
        plt.ylabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.xlabel(r"programmed amplitude*$t_{pulse}$ / $\mu$s")
        psd.gridandtick(plt.gca())