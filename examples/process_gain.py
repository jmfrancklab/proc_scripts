""""
Process Gain of Receiver Chain
==============================

Two files are required for the following example:

File1 contains the analytic signal acquired on the GDS oscilloscope directly
output by the rf source.
for flake
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
from numpy import r_
import pyspecdata as psd
import pylab as plt
from sympy import symbols, Symbol, latex
import sympy as sp
from scipy.interpolate import CubicSpline

attenuator_dB = 40.021  # attenuator between rf source and receiver chain
data_dir = "ODNP_NMR_comp/noise_tests"
file1 = "240123_power_in_analytic.h5"
file2 = "240123_power_out_analytic.h5"
# Make lists to place the power going into the receiver
# chain and out
P_in = []
P_out = []
# Make list to place the frequencies output by the rf
# source.
P_frq_kHz = []
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
for j, nodename in enumerate(all_node_names):
    # second part of nodename contains output frequency in
    # kHz
    rf_frq_kHz = int(nodename.split("_")[1])
    d = psd.nddata_hdf5(
        file1 + "/" + "%s" % nodename,
        directory=psd.getDATADIR(exp_type=data_dir),
    )
    P_frq_kHz.append(rf_frq_kHz)  # Add output frequency to
    #                              list
    # {{{ fit to complex
    A, omega, phi, t = symbols("A omega phi t", real=True)
    d = psd.lmfitdata(d)
    d.functional_form = A * sp.exp(
        1j * 2 * np.pi * omega * t + 1j * phi
    )
    d.set_guess(
        A=dict(value=5e-2, min=1e-4, max=1),
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
    # instantaneous $V_{p} \rightarrow V_{rms}$
    V = d.output("A") / np.sqrt(2)
    V_sq = abs(V) ** 2  # $V_{rms}^{2}$
    P = V_sq / 50  # $\frac{V_{rms}^2}{50 \Omega} = P$
    P_in.append(P)
# {{{ make spline for power in
P_in_data = psd.nddata(P_in, [-1], ["nu"]).setaxis(
    "nu", P_frq_kHz
)
for_plot = P_in_data.C * 1e6  # for plotting in units of μW
P_in_cs = CubicSpline(
    np.array(P_in_data.getaxis("nu")), P_in_data.data
)
# }}}
# {{{ plot P at input of Receiver Chain
with psd.figlist_var() as fl:
    fl.next("Power Input to Receiver Chain")
    fl.plot(for_plot, "ro")
    P_frq_kHz.sort()
    # make axis of finely spaced frequencies
    P_fine = np.linspace(P_frq_kHz[0], P_frq_kHz[-1], 100)
    for_plot_spline = P_in_cs(P_fine) * 1e6  # for plotting
    #                                         in units of μW
    fl.plot(P_fine, for_plot_spline, color="red")
    plt.xlabel(r"$\nu_{rf}$ / kHz")
    plt.ylabel("Power / $\mu$W")
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
    "nu", P_frq_kHz
)
for_plot = P_out_data.C * 1e6  # plot in units of μW
P_out_cs = CubicSpline(
    np.array(P_out_data.getaxis("nu")), P_out_data.data
)
# }}}
# {{{ Plot P at output of Receiver Chain
with psd.figlist_var() as fl:
    fl.next("Power Output by Receiver Chain")
    fl.plot(for_plot, "bo")
    for_plot_spline = P_out_cs(P_fine) * 1e6
    fl.plot(P_fine, for_plot_spline, color="blue")
    plt.xlabel(r"$\nu_{rf}$ / kHz")
    plt.ylabel("Power / $\mu$W")
# }}}
# }}}
# {{{ Divide spline function for P in W
gain_fine = P_out_cs(P_fine) / P_in_cs(P_fine)
gain = P_out_data / P_in_data
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
