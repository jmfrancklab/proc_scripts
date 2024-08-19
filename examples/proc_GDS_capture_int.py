import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(colorcyc_list)

V_atten_ratio = 102.35  # attenutation ratio
loose_HH_width = 1e6
int_slop = 1e-6

with psd.figlist_var() as fl:
    for filename, nodename in [
        (
            "240819_amp0p1_beta_max_pulse_capture.h5",
            "pulse_capture_6",
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
                "units weren't set for the t axis or else I can't read them from the hdf5 file!"
            )
            s.set_units("t", "s")
        s *= V_atten_ratio  # attenutation ratio
        s /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)
        fl.next(r"$\sqrt{P_{analytic}}$ vs $t_{pulse}$")
        raw_color = next(color_cycle)
        fl.plot(s, color=raw_color, label="raw analytic")
        abs_color = next(color_cycle)
        fl.plot(abs(s), color=abs_color, label="abs(analytic)")
        # {{{ apply frequency filter
        SW = 1 / (s["t"][1] - s["t"][0])  # sample rate
        carrier = (
            s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
        )  # signal frequency
        fn = SW / 2  # nyquist frequency
        n = np.floor(
            (carrier + fn) / SW
        )  # how far is the carrier from the left side of the spectrum (which is at SW/2), in integral multiples of SW
        nu_a = (
            carrier - n * SW
        )  # find the aliased peak -- again, measuring from the left side
        center = SW - abs(nu_a)
        s.ft("t")
        fl.next("Frequency Domain")
        fl.plot(s, color=raw_color, label="raw analytic")
        fl.plot(abs(s), color=abs_color, label="abs(analytic)")
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
        # }}}
        s.ft("t")
        wide_HH = s.C
        wide_HH["t" : (0, center - loose_HH_width)] *= 0
        wide_HH["t" : (center + loose_HH_width, None)] *= 0
        # }}}
        # {{{ plot application of all filters
        for filtered_data, label, ax_place in [
            (wide_HH, "loose HH", 0.5),
            (for_L, "lorentzian", -0.5),
        ]:
            color = next(color_cycle)
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
