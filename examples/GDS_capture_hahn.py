import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np

V_atten_ratio = 102.35  # attenutation ratio
with psd.figlist_var() as fl:
    for filename, label, p90_range, p180_range in [
        (
            "240805_beta_22p4us_amp1_GDS_1atten_actual.h5",
            "amplitude = 1, beta = 22.4",
            (0, 18.38e-6),
            (30e-6, None),
        ),
        (
            "240805_beta_22p4us_amp0p2_GDS_1atten_actual.h5",
            "amplitude = 0.2, beta = 22.4",
            (0, 47e-6),
            (46e-6, None),
        ),
        (
            "240805_beta_22p4us_amp0p1_GDS_1atten_actual.h5",
            "amplitude = 0.1, beta = 22.4",
            (56e-6, 122e-6),
            (110e-6, None),
        ),
        (
            "240805_beta_22p4us_amp0p05_GDS_1atten_actual.h5",
            "amplitude = 0.05, beta = 22.4",
            (20e-6, 104e-6),
            (100e-6, None),
        ),
    ]:
        s = psd.nddata_hdf5(
            filename + "/GDS_capture",
            directory=psd.getDATADIR(exp_type="ODNP_NMR_comp/noise_tests"),
        )
        fl.basename = label
        s.set_units("t", "s")
        if "ch" in s.dimlabels:
            s = s["ch", 0]
            # {{{ make analytic
            s.ft("t", shift=True)
            s = s["t":(0, None)]
            s *= 2
            s["t":0] *= 0.5
            # }}}
            s.ift("t")
        s *= V_atten_ratio  # attenutation ratio
        s /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)
        fig, ax_list = plt.subplots(1, 2, figsize=(10, 4))
        fig.suptitle(fl.basename)
        ax_list[1].set_ylabel(r"$\sqrt{P_{pulse}}$")
        ax_list[0].set_ylabel(r"$\sqrt{P_{pulse}}$")
        ax_list[1].set_ylabel(None)
        ax_list[1].set_title("180 Pulse")
        fl.next("raw capture", fig=fig)
        fl.plot(s["t":p90_range].C, alpha=0.2, color="blue", label = None, ax=ax_list[0])
        fl.plot(s["t":p180_range].C, alpha=0.2, color="blue", label = None, ax=ax_list[1])
        fl.plot(abs(s["t":p90_range]), color="orange", label = None, ax=ax_list[0])
        fl.plot(abs(s["t":p180_range].C), color="orange", label = None, ax=ax_list[1])
        # {{{ apply frequency filter
        dt = s["t"][1] - s["t"][0]
        SW = 1 / dt
        carrier = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
        if int(carrier / SW) % 2 == 0:
            center = carrier % SW
        else:
            center = SW - (carrier % SW)
        s.ft("t")
        L_s = s.C
        x = s.getaxis('t')
        # {{{ lorentzian filter
        Lambda = 15.19e6-14.61e6
        s_max = abs(L_s).data.max()
        L = psd.nddata(s_max/(1+1j*2*(x-carrier)*(1/Lambda)),['t']).setaxis('t',x)
        L.set_units('t',L_s.get_units('t'))
        L_s *= L
        L_s /= s_max
        # }}}
        # {{{ heaviside hat functions
        loose_left = center - 10e6
        loose_right = center + 10e6
        tight_left = center - 0.05e6
        tight_right = center + 0.05e6
        small_HH = s.C
        wide_HH = s.C
        small_HH["t":(0, tight_left)] *= 0
        small_HH["t":(tight_right, None)] *= 0
        wide_HH["t":(0, loose_left)] *= 0
        wide_HH["t":(loose_right, None)] *= 0
        # }}}
        for s, label, color, ax_place in [
                (small_HH, "tight HH", 'magenta',2),
                (wide_HH, "loose HH", 'lime', 1),
                (L_s, "lorentzian", 'red',-1)
                ]:
            s.ift("t")
            s /= np.sqrt(2)
            s.ft('t')
            fl.plot(abs(s),color = color, alpha = 0.5, label = label)
            s.ift('t')
            fl.next("raw capture", fig=fig)
            fl.plot(abs(s["t":p90_range]), ax=ax_list[0], color=color, alpha = 0.5, label = label)
            int_range = abs(s["t":p90_range]).contiguous(
                lambda x: x > 0.01 * s["t":p90_range].max()
            )[0]
            int_range -= 1e-6
            int_range += 1e-6
            ninet = abs(s["t":int_range]).integrate("t").data.item() * 1e6
            ax_list[0].set_ylabel(r"$\sqrt{P_{pulse}}$")
            ax_list[0].set_title("90 Pulse")
            ax_list[0].text(
                p90_range[0] * 1e6,
                ax_place,
                r"$t_{90} \sqrt{P_{tx}}_{%s} = %f s \sqrt{W}$" %(label, ninet),
            )
            ax_list[0].axvline(int_range[0] * 1e6, ls=":", color = color, alpha=0.5)
            ax_list[0].axvline(int_range[1] * 1e6, ls=":", color = color, alpha=0.5)
            fl.plot(abs(s["t":p180_range].C), ax=ax_list[1], color=color, alpha = 0.5, label = label)
            int_range = abs(s["t":p180_range]).contiguous(
                lambda x: x > 0.01 * s["t":p180_range].max()
            )[0]
            int_range -= 1e-6
            int_range += 1e-6
            eight = abs(s["t":int_range].C).integrate("t").data.item() * 1e6
            ax_list[1].set_ylabel(None)
            ax_list[1].set_title("180 Pulse")
            ax_list[1].text(
                p180_range[0] * 1e6,
                ax_place,
                r"$t_{180} \sqrt{P_{tx}}_{%s} = %f s \sqrt{W}$" %(label,eight),
            )
            ax_list[1].axvline(int_range[0] * 1e6, ls=":", color = color, alpha=0.5)
            ax_list[1].axvline(int_range[1] * 1e6, ls=":", color = color, alpha=0.5)
