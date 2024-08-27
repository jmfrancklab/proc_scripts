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
from pyspecProcScripts import find_apparent_anal_freq
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(
    colorcyc_list
)  # this can be done more than once to spin up multiple lists

V_atten_ratio = 102.2  # attenutation ratio
skip_plots = 45  # diagnostic -- set this to None, and there will be no plots
slicewidth = 1e6
typical_180 = 40e-6  # typical beta for a 180 -- it's really important to get pulses in this regime correct

# Note: the linear threshold seems to vary from one amplitude to the next
with psd.figlist_var() as fl:
    for filename, nodename, linear_threshold in [
        (
            "240819_test_amp0p05_calib_pulse_calib.h5",
            "pulse_calib_3",
            150e-6,
        ),
        (
            "240819_amp0p1_calib_pulse_calib.h5",
            "pulse_calib_1",
            270e-6,
        ),
        (
            "240819_amp0p2_calib_repeat_pulse_calib.h5",
            "pulse_calib_8",
            310e-6,
        ),
    ]:
        d = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        assert (
            d.get_prop("postproc_type") == "GDS_capture_v1"
        ), "The wrong postproc_type was set so you most likely used the wrong script for acquisition"
        amplitude = d.get_prop("acq_params")["amplitude"]
        fl.basename = f"amplitude = {amplitude:.1f}"
        # {{{ ensure units are set, the preprocessing will through out the statement that they were not set but to still be able to run this script we include the following
        if not d.get_units("t") == "s":
            d.set_units("t", "s")
        if not d.get_units("t_pulse") == "s":
            d.set_units("t_pulse", "s")
        # }}}
        d *= V_atten_ratio
        d /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)
        t_pulses = (
            d.getaxis("t_pulse") / 1e-6
        )  # needed for titles of plots when t_pulse is already indexed

        # {{{ functions that streamline plotting the desired pulse length datasets
        def switch_to_plot(d, thisj):
            thislen = t_pulses[thisj]
            fl.next(f"pulse length = {thislen:.2f} μs")

        def indiv_plots(d, thislabel, thiscolor, thisj=None):
            if skip_plots is None:
                return
            # {{{ this will handle both when t_pulse is (e.g showing integration bounded data) and is not indexed already
            j = (
                range(len(d["t_pulse"]))
                if "t_pulse" in d.dimlabels
                else [thisj]
            )
            for idx in j:
                if idx % skip_plots == 0:
                    switch_to_plot(d, idx)
                    fl.plot(
                        d["t_pulse", idx] if "t_pulse" in d.dimlabels else d,
                        alpha=0.2,
                        color=thiscolor,
                        label=thislabel,
                    )

                    plt.ylabel(r"$\sqrt{P}$ / $\mathrm{\sqrt{W}}$")

        # }}}

        # {{{ data is already analytic, and downsampled to below 24 MHz
        indiv_plots(abs(d), "absolute analytic", "orange")  # sqrt(P_p)
        d, nu_a_MHz, _ = find_apparent_anal_freq(d)  # find frequency of signal
        d.ft("t")
        # {{{ Diagnostic to ensure the frequency was properly identified
        fl.next("Frequency Domain")
        fl.plot(d)
        plt.text(
            x=0.5,
            y=0.5,
            s=rf"$\nu_a={nu_a_MHz/1e6:0.2f}$ MHz",
            transform=plt.gca().transAxes,
        )
        assert (0 > nu_a_MHz + 0.5 * slicewidth) or (
            0 < nu_a_MHz - 0.5 * slicewidth
        ), "unfortunately the region I want to filter includes DC -- this is probably not good, and you should pick a different timescale for your scope so this doesn't happen"
        # }}}
        # {{{ apply HH frequency filter
        d["t" : (0, nu_a_MHz - 0.5 * slicewidth)] *= 0
        d["t" : (nu_a_MHz + 0.5 * slicewidth, None)] *= 0
        # }}}
        d.ift("t")
        indiv_plots(abs(d), "filtered analytic", "red")
        # }}}
        thiscolor = next(color_cycle)
        # {{{ set up shape of beta to drop the calculated values into
        beta = d.shape.pop("t").alloc(dtype=np.float64)
        beta.copy_axes(d).set_units(r"s√W")
        # }}}
        for j in range(len(d["t_pulse"])):
            s = d["t_pulse", j]
            int_range = abs(s).contiguous(lambda x: x > 0.05 * s.max())[0]
            # slightly expand int range to include rising edges
            int_range[0] -= 5e-6
            int_range[-1] += 5e-6
            beta["t_pulse", j] = (
                abs(s["t":int_range]).integrate("t").data.item()
            )
            beta["t_pulse", j] /= np.sqrt(2)  # t*sqrt(Prms)
            indiv_plots(
                abs(s["t":int_range]),
                thiscolor="black",
                thislabel="integrated slice",
                thisj=j,
            )
        # {{{ show what we observe -- how does β vary with the programmed pulse length
        fl.basename = None  # reset so all amplitudes are on same plots
        fl.next(r"Measured $\beta$ vs A * $t_{pulse}$")
        beta.rename("t_pulse", "$A t_{pulse}$")
        beta.name(r"$\beta$")
        beta[
            "$A t_{pulse}$"
        ] *= amplitude  # we only want the x axis to be multiplied by amplitude for plotting only!
        fl.plot(
            (beta.C / 1e-6),
            color=thiscolor,
            label="amplitude = %f" % amplitude,
        )
        beta[
            "$A t_{pulse}$"
        ] /= amplitude  # need to divide back out for the determination of the coefficients below
        psd.gridandtick(plt.gca())
        beta.rename("$A t_{pulse}$", "t_pulse")
        # }}}
        # {{{ Identify captures that don't increase in beta - don't use
        decreasing_idx = np.nonzero(~(np.diff(beta.data) > 0))[0]
        if (
            len(decreasing_idx) > 0
        ):  # beta doesn't always increase with increasing pulse length
            fl.plot(
                beta["t_pulse", : decreasing_idx[-1] + 1],
                "x",
                color="r",
                label="can't use these",
            )
            beta = beta["t_pulse", decreasing_idx[-1] + 1 :]
        # }}}
        plt.axhline(  # the linear threshold is the threshold above which beta is linear
            y=linear_threshold / 1e-6,  # as above
            color=thiscolor,
            label=f"linear threshold for amp={amplitude}",
        )
        # {{{ flip data so beta is on x axis now and A*t_pulse is the y axis
        t_us_v_beta = (
            beta.shape.alloc(dtype=np.float64)
            .rename("t_pulse", r"$\beta$")
            .setaxis(r"$\beta$", beta.data)
        )
        t_us_v_beta.data[:] = (
            beta["t_pulse"].copy() / 1e-6
        )  # because our ppg wants μs
        t_us_v_beta.set_units("μs").set_units(r"$\beta$", "s√W")
        # use as temp for ultimate coeff
        # }}}
        # {{{ Determine linear and nonliear coefficients
        c_nonlinear = t_us_v_beta[r"$\beta$":(None, linear_threshold)].C
        c_nonlinear[
            r"$\beta$"
        ] -= linear_threshold  # Taylor expand around the linear threshold rather than 0
        c_nonlinear = c_nonlinear.polyfit(r"$\beta$", order=10)
        c_linear = t_us_v_beta[r"$\beta$":(linear_threshold, None)].polyfit(
            r"$\beta$", order=1
        )
        print(
            f"\n**************** Coefficients for {amplitude} ****************\n"
        )
        print("Non-linear regime coefficients:\n", c_nonlinear)
        print("Linear regime coefficients:\n", c_linear)

        # }}}
        def prog_plen(desired):
            """function that takes the coefficients of the linear and nonlinear
            regions and applies the fit respectively
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
        t_us_v_beta.name("A * $t_{pulse}$")
        fl.plot(
            t_us_v_beta * amplitude,
            ".",
            alpha=0.5,
            label=f"data for {amplitude}",
        )
        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$, zoomed")
        fl.plot(
            t_us_v_beta[r"$\beta$":(None, typical_180)] * amplitude,
            ".",
            alpha=0.5,
            label=f"data for {amplitude}",
        )
        # {{{ we extrapolate past the edges of the data to show how the
        #     nonlinear is poorly behaved for large beta values
        for_extrap = (
            psd.nddata(
                np.linspace(
                    0.5e-6, t_us_v_beta[r"$\beta$"].max() + 10e-6, 500
                ),
                r"$\beta$",
            )
            .set_units("μs")
            .set_units(r"$\beta$", "s√W")
        )
        for_extrap.copy_props(t_us_v_beta)
        fl.plot(
            for_extrap.eval_poly(c_linear, r"$\beta$")[
                r"$\beta$":(None, typical_180)
            ]
            * amplitude,
            "--",
            alpha=0.25,
            label="linear",
        )
        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$")
        fl.plot(
            for_extrap.eval_poly(c_linear, r"$\beta$") * amplitude,
            "--",
            alpha=0.25,
            label="linear",
        )
        full_fit = for_extrap.fromaxis(r"$\beta$").run(prog_plen)
        fl.plot(full_fit * amplitude, alpha=0.5, label="fit")
        plt.axvline(
            x=linear_threshold / 1e-6,  # units of μs
            alpha=0.1,
            color=for_extrap.get_plot_color(),
        )
        fl.next(r"Amplitude*$t_{pulse}$ vs $\beta$, zoomed")
        fl.plot(
            full_fit[r"$\beta$":(None, typical_180)] * amplitude,
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
