r"""
Calculate actual beta as a function of pulse length
===================================================
Assuming the data is the capture of the pulse sequence as seen on the GDS
oscilloscope (acquired using FLInst/examples/calib_pulses.py), 
here the data is converted to analytic power, frequency filtered
and the absolute is taken prior to integrating to return the beta where
:math:`\beta = \frac{1}{\sqrt{2}} \int \sqrt{P(t)} dt` 
"""
import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(
    colorcyc_list
)  # this can be done more than once to spin up multiple lists

V_atten_ratio = 102.2  # attenutation ratio
skip_plots = 33  # diagnostic -- set this to None, and there will be no plots
slicewidth = 1e6
typical_180 = 40e-6  # typical beta for a 180 -- it's really important to get pulses in this regime correct


with psd.figlist_var() as fl:
    for filename, nodename, linear_threshold in [
        (
            "240819_test_amp0p05_calib_pulse_calib.h5",
            "pulse_calib_3",
            150e-6,
        ),  # low power
        (
            "240819_amp0p1_calib_pulse_calib.h5",
            "pulse_calib_1",
            270e-6,
        ),  # low power
        (
            "240819_amp0p2_calib_repeat_pulse_calib.h5",
            "pulse_calib_8",
            310e-6,
        ),  # low power
    ]:
        d = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        assert (
            d.get_prop("postproc_type") == "GDS_capture_v1"
        ), "The wrong postproc_type was set so you most likely used the wrong script for acquisition"
        amplitude = d.get_prop("acq_params")["amplitude"]
        fl.basename = f"amplitude = {amplitude}"
        if not d.get_units("t") == "s":
            print(
                "units weren't set for the t axis or else I can't read them from the hdf5 file!"
            )
            d.set_units("t", "s")
        d *= V_atten_ratio
        d /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)

        def switch_to_plot(d, j):
            thislen = d["t_pulse"][j]
            fl.next(f"pulse length = {thislen}")

        def indiv_plots(d, thislabel, thiscolor):
            if skip_plots is None:
                return
            for j in range(len(d["t_pulse"])):
                if j % skip_plots == 0:
                    switch_to_plot(d, j)
                    fl.plot(
                        d["t_pulse", j],
                        alpha=0.2,
                        color=thiscolor,
                        label=thislabel,
                    )

        indiv_plots(d, "raw", "blue")
        # {{{ data is already analytic, and downsampled to below 24 MHz
        indiv_plots(abs(d), "analytic", "orange")
        d.ft("t")
        # {{{ apply frequency filter
        d.ift("t")
        SW = 1 / (d["t"][1] - d["t"][0])  # sample rate
        carrier = (
            d.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
        )  # signal frequency
        fn = SW / 2  # nyquist frequency
        n = np.floor(
            (carrier + fn) / SW
        )  # how far is the carrier from the left side of the spectrum (which is at SW/2), in integral multiples of SW
        nu_a = (
            carrier - n * SW
        )  # find the aliased peak -- again, measuring from the left side
        center = SW - abs(nu_a)
        d.ft("t")
        d["t" : (0, center - 0.5 * slicewidth)] *= 0
        d["t" : (center + 0.5 * slicewidth, None)] *= 0
        # }}}
        d.ift("t")
        indiv_plots(abs(d), "filtered analytic", "red")
        fl.next("collect filtered analytic", legend=True)
        for j in range(d.shape["t_pulse"]):
            s = d["t_pulse", j].C
            s["t"] -= abs(s).contiguous(lambda x: x > 0.05 * s.max())[0][0]
            fl.plot(abs(s), alpha=0.3, label=f"{j}")
        # }}}
        thislabel = "amplitude = %f" % amplitude
        thiscolor = next(color_cycle)
        beta = d.shape.pop("t").alloc(dtype=np.float64)
        beta.copy_axes(d)
        beta.set_units(r"s√W").set_units("t_pulse", "s")
        for j in range(len(d["t_pulse"])):
            s = d["t_pulse", j]
            thislen = d["t_pulse"][j]
            int_range = abs(s).contiguous(lambda x: x > 0.05 * s.max())[0]
            # slightly expand int range to include rising edges
            int_range[0] -= 5e-6
            int_range[-1] += 5e-6
            beta["t_pulse", j] = (
                abs(s["t":int_range]).integrate("t").data.item()
            )
            beta["t_pulse", j] /= np.sqrt(2)  # t*sqrt(Prms)
            if skip_plots is not None and j % skip_plots == 0:
                switch_to_plot(d, j)
                fl.plot(
                    abs(s["t":int_range]),
                    color="black",
                    label="integrated slice",
                )
                plt.ylabel(r"$\sqrt{P_{pulse}}$")
                plt.text(
                    int_range[0] * 1e6 - 1,
                    -1,
                    r"$t_{90} \sqrt{P_{tx}} = %f \mathrm{μs} \sqrt{\mathrm{W}}$"
                    % (beta["t_pulse", j].item() / 1e-6),
                )
        # {{{ show what we observe -- how does β vary with the programmed pulse length
        fl.basename = None
        fl.next(r"Measured $\beta$ vs A * $t_{pulse}$")
        beta.rename("t_pulse", "$A t_{pulse}$")
        beta.name(r"$\beta$")
        fl.plot(
            (beta.C / 1e-6).set_units("μs√W") * amplitude,
            color=thiscolor,
            label=thislabel,
        )
        psd.gridandtick(plt.gca())
        beta.rename("$A t_{pulse}$", "t_pulse")
        # }}}
        decreasing_idx = np.nonzero(~(np.diff(beta.data) > 0))[0]
        if (
            len(decreasing_idx) > 0
        ):  # beta doesn't always increase with increasing pulse length
            fl.plot(
                beta["t_pulse", : decreasing_idx[-1] + 1] * amplitude,
                "x",
                color="r",
                label="can't use these",
            )
            beta = beta["t_pulse", decreasing_idx[-1] + 1 :]
        plt.axvline(  # the linear threshold is the threshold above which beta is linear
            x=linear_threshold / 1e-6,  # as above
            color=thiscolor,
            label=f"linear threshold for amp={amplitude}",
        )
        t_us_v_beta = beta.shape.alloc(dtype=np.float64).rename(
            "t_pulse", "beta"
        )
        t_us_v_beta.setaxis("beta", beta.data)
        t_us_v_beta.data[:] = (
            beta["t_pulse"].copy() / 1e-6
        )  # because our ppg wants μs
        t_us_v_beta.set_units("μs").set_units("beta", "s√W")
        # use as temp for ultimate coeff
        c_nonlinear = t_us_v_beta["beta":(None, linear_threshold)].C
        c_nonlinear[
            "beta"
        ] -= linear_threshold  # Taylor expand around the linear threshold rather than 0
        c_nonlinear = c_nonlinear.polyfit("beta", order=10)
        c_linear = t_us_v_beta["beta":(linear_threshold, None)].polyfit(
            "beta", order=1
        )
        print(
            "Non-linear regime coefficients for %s:" % fl.basename, c_nonlinear
        )
        print("Linear regime coefficients for %s:" % fl.basename, c_linear)

        def prog_plen(desired):
            """function that takes the coefficients of the linear and nonlinear
            regions and applies the fit respectively to calculate the pulse time that
            will return the desired beta value
            """

            def zonefit(desired):
                if desired > linear_threshold:
                    return np.polyval(c_linear[::-1], desired)
                else:
                    return np.polyval(
                        c_nonlinear[::-1], desired - linear_threshold
                    )

            ret_val = np.vectorize(zonefit)(desired)
            if ret_val.size > 1:
                return ret_val
            else:
                return ret_val.item()

        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$", legend=True)
        t_us_v_beta.set_plot_color_next()
        fl.plot(
            t_us_v_beta * amplitude,
            ".",
            alpha=0.5,
            label="data for %s" % thislabel,
        )
        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$, zoomed")
        fl.plot(
            t_us_v_beta["beta":(None, typical_180)] * amplitude,
            ".",
            alpha=0.5,
            label="data for %s" % thislabel,
        )
        # {{{ we extrapolate past the edges of the data to show how the
        #     nonlinear is poorly behaved for large beta values
        for_extrap = (
            psd.nddata(
                np.linspace(0.5e-6, t_us_v_beta["beta"].max() + 10e-6, 500),
                "beta",
            )
            .set_units("μs")
            .set_units("beta", "s√W")
        )
        for_extrap.copy_props(t_us_v_beta)
        fl.plot(
            for_extrap.eval_poly(c_linear, "beta")["beta":(None, typical_180)]
            * amplitude,
            "--",
            alpha=0.25,
            label="linear",
        )
        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$")
        fl.plot(
            for_extrap.eval_poly(c_linear, "beta") * amplitude,
            "--",
            alpha=0.25,
            label="linear",
        )
        full_fit = for_extrap.fromaxis("beta").run(prog_plen)
        fl.plot(full_fit * amplitude, alpha=0.5, label="fit")
        plt.axvline(
            x=linear_threshold / 1e-6,  # units of μs
            alpha=0.1,
            color=for_extrap.get_plot_color(),
        )
        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$, zoomed")
        fl.plot(
            full_fit["beta":(None, typical_180)] * amplitude,
            alpha=0.5,
            label="fit",
        )
        for j in (
            r"Amplitude*$t_{pulse}$ vs $\beta$",
            r"Amplitude*$t_{pulse}$ vs $\beta$, zoomed",
        ):
            fl.next(j)
            psd.gridandtick(plt.gca())
            plt.ylabel(r"$At_{pulse}$ / $\mathrm{\mu s}$")
            plt.xlabel(r"$\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        # }}}
