import pyspecdata as psd
import matplotlib.pyplot as plt
from numpy import sqrt

with psd.figlist_var() as fl:
    for filename, label, p90_range, p180_range, p90_int, p180_int in [
        (
            "240722_50tau_p90_11_GDS_1atten.h5",
            "240710 ampllitude = 1 \n p90 = 11",
            (73e-6, 120e-6),
            (127e-6, 205e-6),
            (88e-6, 103.5e-6),
            (152e-6, 177.5e-6),
        ),
    ]:
        s = psd.nddata_hdf5(
            filename + "/GDS_capture",
            directory=psd.getDATADIR(exp_type="ODNP_NMR_comp/noise_tests"),
        )
        fl.basename = label
        s.set_units("t", "s")
        s = s["ch", 0]
        s *= 102  # attenutation ratio
        s /= sqrt(2)  # Vrms
        s /= sqrt(50)  # V/sqrt(R) = sqrt(P)
        fig, ax_list = plt.subplots(1, 2, figsize=(10, 4))
        fig.suptitle(fl.basename)
        # {{{ Plot raw capture
        ax_list[1].set_ylabel(r"$\sqrt{P_{pulse}}$")
        ax_list[0].set_ylabel(r"$\sqrt{P_{pulse}}$")
        ax_list[1].set_ylabel(None)
        ax_list[1].set_title("180 Pulse")
        fl.next("raw capture x 102 for attenuation ratio")
        fl.plot(s["t":p90_range].C, alpha=0.2, color="blue", ax=ax_list[0])
        fl.plot(s["t":p180_range].C, alpha=0.2, color="blue", ax=ax_list[1])
        # }}}
        # {{{ make analytic
        s.ft("t", shift=True)
        s = s["t":(0, None)]
        s *= 2
        s["t":0] *= 0.5
        # }}}
        s.ift("t")
        # {{{ Plot analytic
        fl.plot(abs(s["t":p90_range]), color="orange", ax=ax_list[0])
        fl.plot(abs(s["t":p180_range].C), color="orange", ax=ax_list[1])
        s.ft("t")
        # }}}
        # {{{ frequency filter
        s["t":(0, 11e6)] *= 0
        s["t":(24e6, None)] *= 0
        # }}}
        s.ift("t")
        # {{{ plot and integrate 90 pulse
        fl.plot(abs(s["t":p90_range]), ax=ax_list[0], color="red")
        ninet = abs(s["t":p90_int].C).integrate("t")
        ninet *= 1e6
        ax_list[0].set_ylabel(r"$\sqrt{P_{pulse}}$")
        ax_list[0].set_title("90 Pulse")
        ax_list[0].text(
            p90_range[0] * 1e6,
            -5,
            r"$t_{90} \sqrt{P_{tx}} = %f s \sqrt{W}$" % ninet.data.item(),
        )
        ax_list[0].axvline(p90_int[0] * 1e6, ls=":", alpha=0.2)
        ax_list[0].axvline(p90_int[1] * 1e6, ls=":", alpha=0.2)
        # }}}
        # {{{ plot and integrate 180 pulse
        eight = abs(s["t":p180_int].C).integrate("t")
        eight *= 1e6
        fl.plot(abs(s["t":p180_range].C), ax=ax_list[1], color="red")
        ax_list[1].set_ylabel(None)
        ax_list[1].set_title("180 Pulse")
        ax_list[1].text(
            p180_range[0] * 1e6,
            -5,
            r"$t_{180} \sqrt{P_{tx}} = %f s \sqrt{W}$" % eight.data.item(),
        )
        ax_list[1].axvline(p180_int[0] * 1e6, ls=":", alpha=0.2)
        ax_list[1].axvline(p180_int[1] * 1e6, ls=":", alpha=0.2)
        # }}}
