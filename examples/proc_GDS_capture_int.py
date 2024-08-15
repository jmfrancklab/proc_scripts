import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np

V_atten_ratio = 102.35  # attenutation ratio
loose_HH_width = 5e6
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
        fn = SW / 2
        m = np.round(carrier / SW)
        nu_a = (-(1**m)) * (carrier - 2 * m * fn)
        center = SW - abs(nu_a)
        # }}}
        s.ft("t")
        wide_HH = s.C
        wide_HH["t" : (0, center - loose_HH_width)] *= 0
        wide_HH["t" : (center + loose_HH_width, None)] *= 0
        # }}}
        # {{{ plot application of all filters
        for filtered_data, label, color, ax_place in [
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
            beta = (
                abs(filtered_data).integrate("t").data.item()
                / np.sqrt(2)
                * 1e6
            )
            plt.text(
                1,
                ax_place,
                r"$\beta_{%s} = %f \mu s \sqrt{W}$" % (label, beta),
            )
        plt.ylabel(r"$\sqrt{P}$ / $\mathrm{\sqrt{W}}$")
        # }}}