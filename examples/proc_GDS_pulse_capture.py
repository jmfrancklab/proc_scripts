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
from numpy import sqrt, array

show_all = False
with psd.figlist_var() as fl:
    for filename, nodename, beta90_range, beta180_range in [
        (
            "240729_test_p90_pulse_calib.h5",
            "pulse_calib_2",
            (None, 34e-6),
            (34e-6, None),
        ),
    ]:
        d = psd.nddata_hdf5(
            filename + f"/{nodename}",
            directory=psd.getDATADIR(exp_type="ODNP_NMR_comp/test_equipment"),
        )
        beta_1 = []
        beta_2 = []
        for j in range(len(d["p_90"])):
            s = d["p_90", j].C
            s.set_units("t", "s")
            s *= 101.35  # attenutation ratio
            s /= sqrt(2)  # Vrms
            s /= sqrt(50)  # V/sqrt(R) = sqrt(P)
            if show_all:
                fl.basename = r"programmed $\beta$ = %s" % str(
                    s.get_prop("set_beta")[j]
                )
                fig, ax_list = plt.subplots(1, 2, figsize=(10, 4))
                fig.suptitle(fl.basename)
                # {{{ Plot raw capture
                ax_list[1].set_ylabel(r"$\sqrt{P_{pulse}}$")
                ax_list[0].set_ylabel(r"$\sqrt{P_{pulse}}$")
                ax_list[1].set_ylabel(None)
                ax_list[1].set_title("180 Pulse")
                fl.next("Pulse Calibration")
                fl.plot(s["t":beta90_range], alpha=0.2, color="blue", ax=ax_list[0])
                fl.plot(s["t":beta180_range], alpha=0.2, color="blue", ax=ax_list[1])
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
                fl.plot(abs(s["t":beta90_range]), color="orange", ax=ax_list[0])
                fl.plot(abs(s["t":beta180_range]), color="orange", ax=ax_list[1])
                s.ft("t")
                # }}}
            # {{{ frequency filter
            s["t":(0, 11e6)] *= 0
            s["t":(24e6, None)] *= 0
            # }}}
            s.ift("t")
            if show_all:
                # {{{ plot and integrate 90 pulse
                fl.plot(abs(s["t":beta90_range]), ax=ax_list[0], color="red")
                fl.plot(abs(s["t":beta180_range]), ax=ax_list[1], color="red")
                # }}}
            beta90_int = abs(s["t":beta90_range]).contiguous(
                lambda x: x > 0.01 * s["t":beta90_range].max()
            )[0]
            beta180_int = abs(s["t":beta180_range]).contiguous(
                lambda x: x > 0.01 * s["t":beta180_range].max()
            )[0]
            beta1 = abs(s["t":beta90_int]).integrate("t").data.item() * 1e6
            beta_1.append(beta1)
            beta2 = abs(s["t":beta180_int]).integrate("t").data.item() * 1e6
            beta_2.append(beta2)
            if show_all:
                ax_list[0].set_ylabel(r"$\sqrt{P_{pulse}}$")
                ax_list[0].set_title("90 Pulse")
                ax_list[0].text(
                    beta90_int[0] * 1e6 - 1,
                    -1,
                    r"$t_{90} \sqrt{P_{tx}} = %f s \sqrt{W}$" % beta1,
                )
                ax_list[0].axvline(beta90_int[0] * 1e6, ls=":", alpha=0.2)
                ax_list[0].axvline(beta90_int[1] * 1e6, ls=":", alpha=0.2)
                # }}}
                ax_list[1].set_ylabel(None)
                ax_list[1].set_title("180 Pulse")
                ax_list[1].text(
                    beta180_range[0] * 1e6,
                    -1,
                    r"$t_{180} \sqrt{P_{tx}} = %f s \sqrt{W}$" % beta2,
                )
                ax_list[1].axvline(beta180_int[0] * 1e6, ls=":", alpha=0.2)
                ax_list[1].axvline(beta180_int[1] * 1e6, ls=":", alpha=0.2)
                # }}}
        # {{{ make nddata for beta of 90 pulse and 180 pulse
        if show_all:
            fl.basename = None
        fig2, axs = plt.subplots(1, 2, figsize=(10, 4))
        fig2.suptitle(r"Measured $\beta$ vs $t_{90}$")
        fl.next(r"Measured $\beta$ vs $t_{90}$")
        beta_a = psd.nddata(array(beta_1), ["prog_p90"])
        beta_a.setaxis("prog_p90", d.get_prop("set_p90s"))
        beta_b = psd.nddata(array(beta_2), ["prog_p90"])
        beta_b.setaxis("prog_p90", d.get_prop("set_p180s"))
        fl.plot(beta_a, "o", ax=axs[0])
        axs[0].set_ylabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        axs[0].set_xlabel(r"programmed $t_{90}$ / $\mathrm{\mu s}$")
        axs[0].set_title("First Pulse")
        fl.plot(beta_b, "o", ax=axs[1])
        axs[1].set_ylabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        axs[1].set_xlabel(r"programmed $t_{90}$ / $\mathrm{\mu s}$")
        axs[1].set_title("Second Pulse")
