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

attenuator_dB = 40.021  # exact (measured) attenuation of
#                         attenuation assembly between
#                         source and receiver chain
data_dir = "ODNP_NMR_comp/noise_tests"
file1 = "240123_power_in_analytic.h5"
file2 = "240123_power_out_analytic.h5"
# Nodenames for both files match so extract them here for both files
all_node_names = sorted(
    psd.find_file(
        file1,
        exp_type=data_dir,
        return_list=True,
    ),
    key=lambda x: int(x.split("_")[1]),
)
input_frq_kHz = np.array([int(j.split("_")[1]) for j in all_node_names])
# {{{ Allocate an nddata that's big enough to just
#     store the data for all frequencies
input_data = (
    psd.ndshape([len(all_node_names)], ["nu"]).alloc().set_units("nu", "Hz")
)
input_data["nu"] = input_frq_kHz * 1e3
input_data.rename("nu",r"$\nu$").name("Input Power").set_units("W")
# Copy shape of data but rename for output power
output_data = input_data.C
output_data.name("Output Power").set_units("W")
# }}}
# {{{ set symbols and function for fit
A, omega, phi, t = symbols("A omega phi t", real=True)
fit_function = A * sp.exp(1j * 2 * np.pi * omega * t + 1j * phi)
# }}}
with psd.figlist_var() as fl:
    # {{{ Calculate Power Going into Receiver Chain as Function
    #     of Frequency
    for j, nodename in enumerate(all_node_names):
        # second part of nodename contains output frequency in
        # kHz
        d = psd.nddata_hdf5(
            file1 + "/" + "%s" % nodename,
            directory=psd.getDATADIR(exp_type=data_dir),
        )
        # {{{ fit to complex
        d = psd.lmfitdata(d)
        d.functional_form = fit_function
        d.set_guess(
            A=dict(value=5e-2, min=1e-4, max=1),
            omega=dict(
                value=input_data[r"$\nu$"][j],
                min=input_data[r"$\nu$"][j] - 1e4,
                max=input_data[r"$\nu$"][j] + 1e4,
            ),
            phi=dict(value=0.75, min=-np.pi, max=np.pi),
        )
        d.fit(use_jacobian=False)
        d.eval()
        # }}}
        # calculate (cycle averaged) power from amplitude of the
        # analytic signal:
        input_data[r"$\nu$", j] = abs(d.output("A")) ** 2 / 2 / 50
    input_data.sort(r"$\nu$")
    # {{{ plot P at input of Receiver Chain
    fl.next("Power Input to Receiver Chain")
    input_data.set_plot_color("r")
    input_data.human_units(scale_data=True)
    input_spline = input_data.spline_lambda()
    P_fine = np.linspace(
        input_data[r"$\nu$"][0],
        input_data[r"$\nu$"][-1],
        500,
    )
    fl.plot(input_spline(P_fine))
    fl.plot(input_data, "o")
    # }}}
    # }}}
    # {{{ Calculate Power exiting the Receiver Chain as Function
    #     of Frequency
    for j, nodename in enumerate(all_node_names):
        rf_frq_kHz = int(nodename.split("_")[1])
        d = psd.nddata_hdf5(
            file2 + "/" + "%s" % nodename,
            directory=psd.getDATADIR(exp_type=data_dir),
        )
        # {{{ fit to complex
        d = psd.lmfitdata(d)
        d.functional_form = fit_function
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
        d.eval()
        # }}}
        # Calculate (cycle averages) power from amplitude of the
        # analytic signal
        output_data[r"$\nu$", j] = abs(d.output("A")) ** 2 / 2 / 50
    output_data.sort(r"$\nu$")
    # {{{ Plot P at output of Receiver Chain
    fl.next("Power Output by Receiver Chain")
    output_data.set_plot_color("b")
    output_data.human_units(scale_data=True)
    output_spline = output_data.spline_lambda()
    fl.plot(output_data, "o")
    fl.plot(output_spline(P_fine))
    # }}}
    # }}}
    # {{{ Calculate and plot gain
    gain_fine_dB = (
        10 * np.log10(output_spline(P_fine) / input_spline(P_fine))
        + attenuator_dB
    )
    gain_dB = 10 * np.log10(output_data / input_data) + attenuator_dB
    # {{{ Plot gain
    fl.next("Gain")
    fl.plot(gain_fine_dB, color="purple")
    fl.plot(gain_dB, "o", color="purple")
    plt.ylabel("Gain / dB")
    # }}}
