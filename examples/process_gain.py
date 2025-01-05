""""
Process Gain of Receiver Chain
==============================

Two files are required for the following example:

File1 contains the analytic signal acquired on the GDS oscilloscope directly
output by the rf source.
File2 contains the analytic signal acquired on the GDS oscilloscope at the
output of the receiver chain when the same signal of File1 is injected into a
calibrated attenuator followed by the input of the receiver chain.

Each node pertains to signal with a different frequency (in kHz) so that the
final plots are a function of frequency.
The plots are produced in the following order:
    1) Power output directly from the rf source as a function of frequency
    2) Power output by the receiver chain as a function of frequency
    3) Gain of the receiver chain as a function of frequency
"""

import numpy as np
import pyspecdata as psd
import pylab as plt
from sympy import symbols
import sympy as sp
from scipy.interpolate import CubicSpline

attenuator_dB = 40.021  # exact (measured) attenuation of
#                         attenuation assembly between
#                         source and receiver chain
data_dir = "ODNP_NMR_comp/noise_tests"
file1 = "240123_power_in_analytic.h5"
file2 = "240123_power_out_analytic.h5"
# TODO ☐: you're still using the following, which is sloppy
# -- fix
P_out = []
# {{{ Calculate Power Going into Receiver Chain as Function
#     of Frequency
all_node_names = sorted(
    psd.find_file(
        file1,
        exp_type=data_dir,
        return_list=True,
    ),
    key=lambda x: int(x.split("_")[1]),
)
# TODO ☐: shouldn't the previous and the following be the
#         same for both input and output? Because you know
#         that's how the data is acquired? So, could you
#         save effort by just doing both of these for the
#         one file?
# TODO ☐: In the allocation block, you could then just create output_data as a
#         copy of input_data, as wel
input_frq_kHz = np.array(
    [int(j.split("_")[1]) for j in all_node_names]
)
# {{{ allocate an nddata that's bit enough to just
#     store the data for all frequencies
input_data = (
    psd.ndshape([len(all_node_names)], ["nu"])
    .alloc()
    .set_units("nu", "Hz")
)
input_data["nu"] = input_frq_kHz * 1e3
input_data.name("input power").set_units("W")
# }}}
for j, nodename in enumerate(all_node_names):
    # second part of nodename contains output frequency in
    # kHz
    d = psd.nddata_hdf5(
        file1 + "/" + "%s" % nodename,
        directory=psd.getDATADIR(exp_type=data_dir),
    )
    # {{{ fit to complex
    A, omega, phi, t = symbols("A omega phi t", real=True)
    d = psd.lmfitdata(d)
    d.functional_form = A * sp.exp(
        1j * 2 * np.pi * omega * t + 1j * phi
    )
    d.set_guess(
        A=dict(value=5e-2, min=1e-4, max=1),
        omega=dict(
            value=input_data["nu"][j],
            min=input_data["nu"][j] - 1e4,
            max=input_data["nu"][j] + 1e4,
        ),
        phi=dict(value=0.75, min=-np.pi, max=np.pi),
    )
    d.fit(use_jacobian=False)
    # }}}
    # calculate (cycle averaged) power from amplitude of the
    # analytic signal:
    input_data["nu", j] = abs(d.output("A")) ** 2 / 2 / 50
input_data.sort("nu")
# {{{ plot P at input of Receiver Chain
with psd.figlist_var() as fl:
    fl.next("Power Input to Receiver Chain")
    input_data.set_plot_color("r")
    input_data.human_units(scale_data=True)
    # TODO ☐: be sure you pull most recent pyspecdata -- I
    #         modified spline to copy props so the color
    #         will match
    myspline = input_data.spline_lambda()
    P_fine = np.linspace(
        input_data["nu"][0],
        input_data["nu"][-1],
        500,
    )
    fl.plot(myspline(P_fine))
    fl.plot(input_data, "o")
# }}}
# }}}
# {{{ Calculate Power exiting the Receiver Chain as Function
#     of Frequency
all_node_names = sorted(
    psd.find_file(
        file2,
        exp_type="ODNP_NMR_comp/noise_tests",
        return_list=True,
    ),
    key=lambda x: int(x.split("_")[1]),
)
for j, nodename in enumerate(all_node_names):
    rf_frq_kHz = int(nodename.split("_")[1])
    d = psd.nddata_hdf5(
        file2 + "/" + "%s" % nodename,
        directory=psd.getDATADIR(exp_type=data_dir),
    )
    # {{{ fit to complex
    A, omega, phi, t = symbols("A omega phi t", real=True)
    func_form = A * sp.exp(
        1j * 2 * np.pi * omega * t + 1j * phi
    )
    d = psd.lmfitdata(d)
    d.functional_form = func_form
    d.set_guess(
        A=dict(value=15e-2, min=1e-4, max=1),
        omega=dict(
            value=rf_frq_kHz * 1e3,
            min=rf_frq_kHz * 1e3 - 1e4,
            max=rf_frq_kHz * 1e3 + 1e4,
        ),
        phi=dict(value=0.75, min=-np.pi, max=np.pi),
    )
    d.fit(use_jacobian=False)
    fit = d.eval()
    # }}}
    # {{{ convert to Power
    V = d.output("A")
    Vrms = V / np.sqrt(
        2
    )  # instantaneous $V_{p} \rightarrow V_{rms}$
    Vrms_sq = abs(Vrms) ** 2  # $V_{rms}^{2}$
    P = Vrms_sq / 50  # $\frac{V_{rms}^2}{50 \Omega} = P$
    P_out.append(P)
    # }}}
# {{{ make spline for power out
P_out_data = psd.nddata(P_out, [-1], ["nu"]).setaxis(
    "nu", input_frq_kHz
)
for_plot = P_out_data.C * 1e6  # plot in units of μW
P_out_cs = CubicSpline(
    np.array(P_out_data.getaxis("nu")), P_out_data.data
)
# }}}
# {{{ Plot P at output of Receiver Chain
# TODO ☐: below, you have a second (and third!!!) figurelist
#         -- a pretty fundamental nono that we've talked
#         about many times
with psd.figlist_var() as fl:
    fl.next("Power Output by Receiver Chain")
    fl.plot(for_plot, "bo")
    for_plot_spline = P_out_cs(P_fine) * 1e6
    fl.plot(P_fine, for_plot_spline, color="blue")
    plt.xlabel(r"$\nu_{rf}$ / kHz")
    plt.ylabel(r"Power / $\mu$W")
# }}}
# }}}
# {{{ Divide spline function for P in W
gain_fine = P_out_cs(P_fine) / P_in_cs(P_fine)
gain = P_out_data / input_data
# }}}
with psd.figlist_var() as fl:
    fl.next("Gain")
    gain_fine_dB = 10 * np.log10(gain_fine) + attenuator_dB
    fl.plot(P_fine, gain_fine_dB, color="purple")
    # convert gain data points to dB
    gain.data = 10 * np.log10(gain.data) + attenuator_dB
    fl.plot(gain, "o", color="purple")
    plt.xlabel(r"$\nu_{rf}$ / kHz")
    plt.ylabel("Gain / dB")
