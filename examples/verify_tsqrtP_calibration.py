r"""
Verify the pulse calibration
============================
Once we have the pulse-length calibration in place, we want to simply verify that we get the
β values that we program, where
:math:`\beta = \int \sqrt{P(t)} dt` 

This is performed with amplitude pulses driving into an attenuator and measured by a scope,
as for the calibration.
"""
import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(
    colorcyc_list
)  # this can be done more than once to spin up multiple lists

V_atten_ratio = 102.35  # attenutation ratio
skip_plots = 33  # diagnostic -- set this to None, and there will be no plots
linear_threshold = 100e-6
slicewidth = 2e6


with psd.figlist_var() as fl:
    for filename, nodename in [
        (
            "240805_test_calib_amp0p05_pulse_calib.h5",
            "pulse_calib_1",
        ),  # low power
        (
            "240805_test_calib_amp0p1_a_pulse_calib.h5",
            "pulse_calib_2",
        ),  # low power
        (
            "240805_test_calib_amp0p2_a_pulse_calib.h5",
            "pulse_calib_1",
        ),  # low power
        (
            "240805_test_calib_amp1_pulse_calib.h5",
            "pulse_calib_1",
        ),  # low power
    ]:
        d = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        assert (
            d.get_prop("postproc_type") is not None
        ), "No postproc type was set upon acquisition"
        amplitude = d.get_prop("acq_params")["amplitude"]
        fl.basename = f"amplitude = {amplitude}"
        if not d.get_units("t") == "s":
            print(
                "************ AG still needs to finish pyspecdata PR to save units!!! ************"
            )
            d.set_units("t", "s")  # why isn't this done already??
        d *= V_atten_ratio
        d /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)
        pulse_lengths = d.get_prop("programmed_t_pulse_us")

        def switch_to_plot(d, j):
            thislen = pulse_lengths[j]
            fl.next(f"pulse length = {thislen}")

        def indiv_plots(d, thislabel, thiscolor):
            if skip_plots is None:
                return
            for j in range(len(d["beta"])):
                if j % skip_plots == 0:
                    switch_to_plot(d, j)
                    fl.plot(
                        d["beta", j],
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
        n = np.floor(
            (carrier + SW / 2) / SW
        )  # how far is the carrier from the left side of the spectrum (which is at SW/2), in integral multiples of SW
        nu_a = (
            -SW / 2 + (carrier + SW / 2) - n * SW
        )  # find the aliased peak -- again, measuring from the left side
        center = SW - abs(nu_a)
        d.ft("t")
        d["t" : (0, center - 0.5 * slicewidth)] *= 0
        d["t" : (center + 0.5 * slicewidth, None)] *= 0

        # }}}
        d.ift("t")
        indiv_plots(abs(d), "filtered analytic", "red")
        fl.next("collect filtered analytic", legend=True)
        for j in range(d.shape["beta"]):
            s = d["beta", j].C
            s["t"] -= abs(s).contiguous(lambda x: x > 0.03 * s.max())[0][0]
            fl.plot(abs(s), alpha=0.3, label=f"{j}")
        # }}}
        thislabel = "amplitude = %f" % amplitude
        thiscolor = next(color_cycle)
        beta = d.shape.pop("t").alloc(dtype=np.float64)
        beta.copy_axes(d)
        beta.set_units(r"s√W").set_units("beta", "s√W")
        for j in range(len(d["beta"])):
            s = d["beta", j]
            thislen = d["beta"][j]
            int_range = abs(s).contiguous(lambda x: x > 0.03 * s.max())[0]
            # slightly expand int range to include rising edges
            int_range[0] -= 2e-6
            int_range[-1] += 2e-6
            beta["beta", j] = abs(s["t":int_range]).integrate("t").data.item()
            beta["beta", j] /= np.sqrt(2)  # Vrms
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
                    % (beta["beta", j].item() / 1e-6),
                )
        # {{{ show what we observe -- how does β vary with the programmed pulse length
        fl.basename = None
        fl.next(r"Measured $\beta$ vs programmed $\beta$")
        fl.plot(
            (beta.C / 1e-6).set_units("μs√W"),
            color=thiscolor,
            label=thislabel,
        )
        plt.xlabel(r"Programmed $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.ylabel(r"Measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        psd.gridandtick(plt.gca())
        # }}}
