import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(
    colorcyc_list
)  # this can be done more than once to spin up multiple lists

V_atten_ratio = 102.35  # attenutation ratio
loose_HH_width = 5
tight_HH_width = 1
skip_plots = 1

with psd.figlist_var() as fl:
    for filename, nodename in [
        (
            "240805_test_calib_amp1_pulse_calib.h5",
            "pulse_calib_1",
        ),
    ]:
        s = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        amplitude = s.get_prop("acq_params")["amplitude"]
        beta_us_sqrt_W = s.get_prop("acq_params")["beta_90_s_sqrtW"]
        fl.basename = f"amplitude = {amplitude}, beta = {beta_us_sqrt_W}"
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
        if "beta" in s.dimlabels:
            s = s["beta", -1]
            s.rename("t", "t_pulse")
        s *= V_atten_ratio  # attenutation ratio
        s /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)

        def switch_to_plot(d, j):
            thislen = s["t_pulse"][j]
            fl.next(f"pulse length = {thislen}")

        def indiv_plots(d, thislabel, thiscolor):
            if skip_plots is None:
                return
            for j in range(len(s["t_pulse"])):
                if j % skip_plots == 0:
                    switch_to_plot(s, j)
                    fl.plot(
                        s["t_pulse", j],
                        alpha=0.2,
                        color=thiscolor,
                        label=thislabel,
                    )

        indiv_plots(s, "raw", "tab:blue")
        indiv_plots(abs(s), "analytic", "tab:orange")
        # {{{ apply frequency filter
        carrier = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
        s.ft("t")
        fl.next("f domain %s" % fl.basename)
        indiv_plots(abs(s), "analytic", "tab:orange")
        x = s.getaxis("t")
        for_L = s.C  # make copy for lorentzian filter
        # {{{ lorentzian filter
        Lambda = 15.19e6 - 14.61e6
        s_max = abs(for_L).data.max()
        L = (
            psd.nddata(s_max / (1 + 1j * 2 * (x - carrier) * (1 / Lambda)), ["t"])
            .setaxis("t", x)
            .set_units("t", for_L.get_units("t"))
        )
        for_L *= L
        for_L /= s_max  # for plotting
        # }}}
        # {{{ heaviside hat functions
        # {{{ determine center frequency
        s.ift("t")
        dt = s["t"][1] - s["t"][0]
        SW = 1 / dt
        if int(carrier / SW) % 2 == 0:
            center = carrier % SW
        else:
            center = SW - (carrier % SW)
        # }}}
        center /= 1e6
        print(center)
        s.ft("t")
        small_HH = s.C
        wide_HH = s.C
        small_HH["t" : (0, center - tight_HH_width)] *= 0
        small_HH["t" : (center + tight_HH_width, None)] *= 0
        wide_HH["t" : (0, center - loose_HH_width)] *= 0
        wide_HH["t" : (center + loose_HH_width, None)] *= 0
        # }}}
        # {{{ plot application of all filters
        for filtered_data, label, color, ax_place in [
            (small_HH, "tight HH", "magenta", 2),
            (wide_HH, "loose HH", "lime", 1),
            (for_L, "lorentzian", "red", -1),
        ]:
            filtered_data.ift("t")
            indiv_plots(abs(filtered_data), f"filtered analytic - {label}", f"{color}")
            indiv_plots(abs(s), f"analytic - {label}", "red")
            filtered_data /= np.sqrt(2)
            filtered_data.ft("t")
            fl.next("f domain %s" % filename)
            fl.plot(abs(filtered_data), color=color, alpha=0.5, label=label)
            plt.axvline(center)
            plt.axvline(center - tight_HH_width)
            print(center - tight_HH_width)
            fl.show()
            quit()
            filtered_data.ift("t")
            int_range = abs(filtered_data).contiguous(
                lambda x: x > 0.01 * filtered_data.max()
            )[0]
            int_range -= 1e-6
            int_range += 1e-6
            beta = abs(filtered_data["t":int_range]).integrate("t").data.item() * 1e6
            plt.text(
                int_range[0] * 1e6,
                ax_place,
                r"$t_{90} \sqrt{P_{tx}}_{%s} = %f s \sqrt{W}$" % (label, beta),
            )
            plt.axvline(int_range[0] * 1e6, ls=":", color=color, alpha=0.5)
            plt.axvline(int_range[1] * 1e6, ls=":", color=color, alpha=0.5)
