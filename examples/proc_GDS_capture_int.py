import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np

V_atten_ratio = 102.35  # attenutation ratio
loose_HH_width = 5e6
tight_HH_width = 1e6
int_slop = 1e-6

with psd.figlist_var() as fl:
    for filename, nodename in [
        (
            "240812_amp0p1_deblank_50_pulse_capture.h5",
            "pulse_capture_1",
        ),
    ]:
        s = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        amplitude = s.get_prop("acq_params")["amplitude"]
        beta_us_sqrt_W = s.get_prop("acq_params")["beta_90_s_sqrtW"]
        fl.basename = (
            f"amplitude = {amplitude}, $\\beta$ = {beta_us_sqrt_W} \n"
        )
        if not s.get_units("t") == "s":
            print(
                "******** AG still needs to finish pyspecdata PR to save units!!! ********"
            )
            s.set_units("t", "s")
        s *= V_atten_ratio  # attenutation ratio
        s /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)
        fl.next(r"$\sqrt{P_{analytic}}$ vs $t_{pulse}$")
        fl.plot(s, color="tab:blue", label="raw analytic")
        fl.plot(abs(s), color="tab:orange", label="abs(analytic)")
        # {{{ apply frequency filter
        carrier = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
        s.ft("t")
        fl.next("Frequency Domain")
        fl.plot(s, color="tab:blue", label="raw analytic")
        fl.plot(abs(s), color="tab:orange", label="abs(analytic)")
        # {{{ lorentzian filter
        x = s.C.getaxis("t")
        for_L = s.C  # make copy for lorentzian filter
        Lambda = 15.19e6 - 14.61e6
        s_max = abs(for_L).data.max()
        L = (
            psd.nddata(
                s_max / (1 + 1j * 2 * (x - carrier) * (1 / Lambda)), ["t"]
            )
            .setaxis("t", x)
            .set_units("t", for_L.get_units("t"))
        )
        for_L *= L
        for_L /= s_max
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
            (small_HH, "tight HH", "magenta", 1),
            (wide_HH, "loose HH", "lime", 0.5),
            (for_L, "lorentzian", "red", -0.5),
        ]:
            filtered_data.ift("t")
            fl.next(r"$\sqrt{P_{analytic}}$ vs $t_{pulse}$")
            fl.plot(
                abs(filtered_data), color=f"{color}", label=label + " filter"
            )
            filtered_data.ft("t")
            fl.next("Frequency Domain")
            fl.plot(abs(filtered_data), color=color, alpha=0.5, label=label)
            filtered_data.ift("t")
            fl.next(r"$\sqrt{P_{analytic}}$ vs $t_{pulse}$")
            int_range = abs(filtered_data).contiguous(
                lambda x: x > 0.01 * filtered_data.max()
            )[0]
            int_range -= int_slop
            int_range += int_slop
            beta = (
                abs(filtered_data["t":int_range]).integrate("t").data.item()
                / np.sqrt(2)
                * 1e6
            )
            plt.text(
                1,
                ax_place,
                r"$t_{90} \sqrt{P_{tx}}_{%s} = %f \mu s \sqrt{W}$"
                % (label, beta),
            )
            plt.axvline(int_range[0] * 1e6, ls=":", color=color, alpha=0.5)
            plt.axvline(int_range[1] * 1e6, ls=":", color=color, alpha=0.5)
        plt.ylabel(r"$\sqrt{P}$ / $\mathrm{\sqrt{W}}$")
        # }}}
