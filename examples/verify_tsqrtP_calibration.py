r"""
Verify the pulse calibration
============================
Once we have the pulse-length calibration in place, we want to simply verify that we get the
β values that we program, where
:math:`\beta = \frac{1}{\sqrt{2} \int \sqrt{P(t)} dt` 

This is performed with amplitude pulses driving into an attenuator and measured by a scope,
as for the calibration.
"""
import pyspecdata as psd
from pyspecProcScripts import find_apparent_anal_freq
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(
    colorcyc_list
)  # this can be done more than once to spin up multiple lists

V_atten_ratio = 102.2  # attenutation ratio
skip_plots = 33  # diagnostic -- set this to None, and there will be no plots
HH_width = 2e6


with psd.figlist_var() as fl:
    for filename, nodename in [
        (
            "240819_test_amp0p05_calib_pulse_calib.h5",
            "pulse_calib_5",
        ),
        (
            "240819_amp0p1_calib_pulse_calib.h5",
            "pulse_calib_2",
        ),
        (
            "240819_amp0p2_calib_repeat_pulse_calib.h5",
            "pulse_calib_9",
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
        if not d.get_units("t") == "s":
            psd.logger.info(
                psd.strm(
                    "units weren't set for the t axis or else I can't read them from the hdf5 file!"
                )
            )
            d.set_units("t", "s")
        d *= V_atten_ratio  # V at output of amplifier
        d /= np.sqrt(50)  # V/sqrt(R) = sqrt(P)

        # {{{ functions that streamline plotting the desired pulse length datasets
        def switch_to_plot(d, thisj):
            thislen = d.get_prop("programmed_t_pulse")[thisj] / 1e-6
            fl.next(f"pulse length = {thislen:.2f} μs")

        def indiv_plots(d, thislabel, thiscolor, thisj=None):
            if skip_plots is None:
                return
            # {{{ this will handle both when beta is and is not indexed already
            j = range(len(d["beta"])) if "beta" in d.dimlabels else [thisj]
            for idx in j:
                if idx % skip_plots == 0:
                    switch_to_plot(d, idx)
                    fl.plot(
                        d if "beta" not in d.dimlabels else d["beta", idx],
                        alpha=0.2,
                        color=thiscolor,
                        label=thislabel,
                    )
                    plt.ylabel(r"$\sqrt{P}$ / $\mathrm{\sqrt{W}}$")

        # }}}

        # {{{ data is already analytic, and downsampled to below 24 MHz
        indiv_plots(abs(d), "abs(analytic)", "orange")  # sqrt(P_p)
        d, nu_a, _ = find_apparent_anal_freq(d)  # find frequency of signal
        d.ft("t")
        # {{{ Diagnostic to ensure the frequency was properly identified
        fl.next("Frequency Domain")
        fl.plot(d)
        plt.text(
            x=0.5,
            y=0.5,
            s=rf"$\nu_a={nu_a/1e6:0.2f}$ MHz",
            transform=plt.gca().transAxes,
        )
        assert (0 > nu_a * 0.5 * HH_width) or (
            0 < nu_a - 0.5 * HH_width
        ), "unfortunately the region I want to filter includes DC -- this is probablye not good, and you should pick a different timescale for your scope so this doesn't happen"
        # }}}
        # {{{ apply frequency filter
        d["t" : (None, nu_a - 0.5 * HH_width)] *= 0
        d["t" : (nu_a + 0.5 * HH_width, None)] *= 0
        # }}}
        d.ift("t")
        indiv_plots(abs(d), "filtered analytic", "red")
        # }}}
        thiscolor = next(color_cycle)
        # {{{ set up shape of data to drop the correct values in
        verify_beta = d.shape.pop("t").alloc(dtype=np.float64)
        verify_beta.copy_axes(d)
        verify_beta.set_units(r"s√W").set_units("beta", r"s√W")
        # }}}
        for j in range(len(d["beta"])):
            s = d["beta", j]
            int_range = abs(s).contiguous(lambda x: x > 0.03 * s.max())[0]
            # slightly expand int range to include rising edges
            int_range[0] -= 2e-6
            int_range[-1] += 2e-6
            verify_beta["beta", j] = (
                abs(s["t":int_range]).integrate("t").data.item()
            )  # tp * sqrt(P_p)
            verify_beta["beta", j] /= np.sqrt(2)  # tp*sqrt(Prms)
            indiv_plots(
                abs(s["t":int_range]),
                thiscolor="black",
                thislabel="integrated slice",
                thisj=j,
            )
        # {{{ show what we observe -- how does β vary with the programmed pulse length
        fl.basename = None  # we want to plot all amplitudes together now
        fl.next(r"Measured $\beta$ vs programmed $\beta$")
        fl.plot(
            (verify_beta / 1e-6).set_units("μs√W"),
            color=thiscolor,
            label="Amplitude = %0.1f" % amplitude,
        )
        plt.xlabel(r"Programmed $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.ylabel(r"Measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        psd.gridandtick(plt.gca())
        # }}}
