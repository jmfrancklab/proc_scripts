r"""
Calculate actual beta as a function of pulse length
===================================================
Assuming the data is the capture of the pulse sequence as seen on the GDS
oscilloscope (acquired using FLInst/examples/calib_pulses.py), 
here the data is converted to analytic power, frequency filtered
and the absolute is taken prior to integrating to return the beta where
:math:`\beta = \int \sqrt{P(t)} dt` 
"""
import pyspecdata as psd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = cycle(
    colorcyc_list
)  # this can be done more than once to spin up multiple lists

atten_ratio = 101.52  # attenutation ratio
skip_plots = 1  # diagnostic -- set this to None, and there will be no plots
with psd.figlist_var() as fl:
    for filename, nodename, amplitude in [
        ("240731_amp1_calib_fin_pulse_calib.h5", "pulse_calib_1", 1.0),
        ("240731_amp0p1_calib_fin_pulse_calib.h5", "pulse_calib_1", 0.1),
    ]:
        fl.basename = f"amplitude = {amplitude}"
        d = psd.find_file(
            filename, expno=nodename, exp_type="ODNP_NMR_comp/test_equipment"
        )
        d["t_pulse"] = np.float64(
                d["t_pulse"])
        d.set_units("t", "s")  # why isn't this done already??
        d *= atten_ratio
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
        # {{{ make analytic, then filter
        d.ft("t", shift=True)
        d = d["t":(0, None)]
        d.ift("t")
        d *= 2
        d["t":0] *= 0.5
        indiv_plots(abs(d), "analytic", "orange")
        d.ft("t")
        d["t":(0, 9.5e6)] *= 0
        d["t":(11.5e6, None)] *= 0
        d.ift("t")
        indiv_plots(abs(d), "filtered analytic", "red")
        fl.next("collect filtered analytic")
        for j in range(d.shape["t_pulse"]):
            s = d["t_pulse", j].C
            s["t"] -= abs(s).contiguous(lambda x: x > 0.01 * s.max())[0][0]
            fl.plot(abs(s), alpha=0.3)
        # }}}
        thislabel = "amplitude = %f" % amplitude
        thiscolor = next(color_cycle)
        beta = d.shape.pop("t").alloc(dtype=np.float64)
        beta.copy_axes(d)
        beta.set_units(r"μs√W")
        for j in range(len(d["t_pulse"])):
            s = d["t_pulse", j]
            thislen = d["t_pulse"][j]
            int_range = abs(s).contiguous(lambda x: x > 0.1 * s.max())[0]
            # slightly expand int range to include rising edges
            int_range[0] -= 1e-6
            int_range[-1] += 1e-6
            beta["t_pulse", j] = (
                abs(s["t":int_range]).integrate("t").data.item() * 1e6
            )
            beta["t_pulse", j] /= np.sqrt(2)  # Vrms
            # PR COMMENT: JF only read to here -- a bunch of stuff above were comments that weren't incorporated or obvious clean code stuff.  Please review from here to the end again
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
                    r"$t_{90} \sqrt{P_{tx}} = %f s \sqrt{W}$"
                    % beta["t_pulse", j].item(),
                )
        # {{{ make nddata for beta of 90 pulse and 180 pulse
        # {{{ beta vs t
        fl.basename = None
        fl.next(r"Measured $\beta$ vs A * $t_{pulse}$")
        beta["t_pulse"] *= amplitude
        beta.rename("t_pulse", "$A t_{pulse}$")
        fl.plot(beta, "o", color=thiscolor, label=thislabel)
        beta.rename("$A t_{pulse}$", "t_pulse")
        beta["t_pulse"] /= amplitude
        t_v_beta = beta.shape.alloc(dtype=np.float64).rename("t_pulse", "beta")
        t_v_beta.setaxis("beta", beta.data)
        t_v_beta.data[:] = beta["t_pulse"].copy()
        if amplitude > 1:
            linear_regime = (7, None)
        else:
            linear_regime = (2, None)
        c_nonlinear = t_v_beta.polyfit("beta", order=10)
        print(c_nonlinear)
        fl.next(r"$t_{pulse}$ vs $\beta$")
        fl.plot(t_v_beta, "o", color = thiscolor)
        t_v_beta.eval_poly(c_nonlinear, "beta")
        fl.plot(t_v_beta, color=thiscolor, ls=":", alpha=0.5)
        psd.gridandtick(plt.gca())
        plt.ylabel(r"measured $\beta$ / $\mathrm{\mu s \sqrt{W}}$")
        plt.xlabel(r"Amplitude*$t_{pulse}$ / $\mu$s")
        # }}}
