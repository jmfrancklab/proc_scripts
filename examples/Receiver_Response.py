""" Calculate and Fit Receiver Response to Sinc
===============================================

Two files are required for the following example:

    File1 contains the analytic signal acquired on the GDS oscilloscope
    directly output from the AFG output. Each node pertains to signal with a
    different frequency (in kHz) which are fit to a complex function to extract
    the power (in :math:`V^{2}`).
    File2 contains the quadrature signal acquired on the receiver when the same
    signal of File1 is injected into it. Each node pertains to signal with a
    different frequency (in kHz) which are converted to a PSD. The power (in
    :math:`dg^{2}`) is calculated from the peak of the convolved PSD.

    The receiver response is then the ratio of :math:`dg^{2}(\\nu)` to
    :math:`V^{2}(\\nu)` which is fit to a sinc function. The final plot shows
    the square of the data with the fit.
"""
import numpy as np
import pyspecdata as psd
from sympy import symbols
import sympy as sp
import re


def for_plot(thisdata):
    """Converts :math:`dg^{2}/V^{2}` to :math:`dg/\\mu V`
    and sets labels accordingly for pretty plotting"""
    # Take sqrt for units of $dg/V$
    thisdata.data = np.sqrt(thisdata.data)
    # Units of $dg/\mu V$
    thisdata /= 1e6
    thisdata.rename("nu_test", r"$\Delta \nu$")
    thisdata.name(r"$\mathrm{dg_{%s\ \mathrm{kHz}}}/ \mu V$" % SW)
    return thisdata


file1 = "240123_10mV_AFG_GDS_5mV_100MSPS_analytic.h5"
file2 = "240117_afg_sc_10mV_3p9kHz_zoom.h5"
lambda_G = 0.4e3  # Width for Gaussian convolution
# There are less nodes acquired for the control case (since we assume it's
# relatively flat) so I need to separately define them
file1_nodes = psd.find_file(
    re.escape(file1),
    exp_type="ODNP_NMR_comp/noise_tests",
    return_list=True,
)
control_frqs_kHz = np.array([float(j.split("_")[1]) for j in file1_nodes])
file2_nodes = psd.find_file(
    re.escape(file2),
    exp_type="ODNP_NMR_comp/noise_tests",
    return_list=True,
)
# $\nu_{test}$
resp_frq_kHz = np.array([float(j.split("_")[1]) for j in file2_nodes])
# {{{ Calculate input power (acquired on Oscilloscope)
#     To determine the input power, just take a capture of the signal in the
#     time domain and fit to a complex exponential
# {{{ Make empty nddata to drop the calculated $V^{2}$ into with corresponding
#     frequency ($\nu$) output by AFG source, and set the frequencies based on
#     the node names
control = psd.ndshape([len(control_frqs_kHz)], ["nu_test"]).alloc()
control.setaxis("nu_test", control_frqs_kHz).set_units("nu_test", "Hz")
# }}}
for j, nodename in enumerate(file1_nodes):
    rf_frq = control["nu_test"][j]
    d = psd.find_file(
        file1,
        expno=nodename,
        exp_type="ODNP_NMR_comp/noise_tests",
    )
    # {{{ fit signal in t domain to complex exponential
    A, omega, phi, t = symbols("A omega phi t", real=True)
    f = psd.lmfitdata(d)
    f.functional_form = A * sp.exp(1j * 2 * np.pi * omega * t + 1j * phi)
    f.set_guess(
        A=dict(value=5e-2, min=1e-4, max=1),
        omega=dict(
            value=rf_frq,
            min=rf_frq - 1e4,
            max=rf_frq + 1e4,
        ),
        phi=dict(value=0.75, min=-np.pi, max=np.pi),
    )
    f.fit(use_jacobian=False)
    fit = f.eval()
    # }}}
    V_amp = f.output("A")
    control["nu_test", j] = V_amp**2
# {{{ make spline for power going into RX box
control.sort("nu_test")
Pin_spline = control.spline_lambda()
# }}}
# }}}
# {{{ Calculate $dg^2$
for j, nodename in enumerate(file2_nodes):
    d = psd.find_file(
        file2,
        exp_type="ODNP_NMR_comp/noise_tests",
        expno=nodename,
    )
    carrier = (
        d.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
    )  # $\nu_{RX,LO}$
    if j == 0:
        # Allocate and array that's shaped like the data acquired
        # for one frequency but add axis to store the frequency of
        # the test signal frequency
        rec_data1 = (d.shape + ("nu_test", len(resp_frq_kHz))).alloc()
        rec_data1.setaxis("t", d.getaxis("t")).set_units("t", "s")
        rec_data1.setaxis("nu_test", resp_frq_kHz).set_units("nu_test", "Hz")
    rec_data1["nu_test", j] = d
    SW = str(d.get_prop("acq_params")["SW_kHz"])
rec_data1.sort("nu_test")
rec_data1.rename("nScans", "capture")  # To be more consistent with the
#                                         oscilloscope data rename the nScans
#                                         dimension
acq_time = np.diff(rec_data1.getaxis("t")[np.r_[0, -1]]).item()
# {{{ Calculate PSD for each frequency
rec_data1.ft("t", shift=True)  # $dg\sqrt{s/Hz}$
rec_data1 = abs(rec_data1) ** 2  # $dg^{2}*s/Hz$
rec_data1.mean("capture")
rec_data1 /= acq_time  # $dg^{2}/Hz$
# Convolve using $\lambda_{G}$ specified above
rec_data1.convolve("t", lambda_G, enforce_causality=False)
rec_data1.ift("t")
rec_data1.run(np.conj)  # Empirically needed to give offset
#                        that increases with field
rec_data1.ft("t")
rec_data1["t"] += carrier
# }}}
# {{{ Calculate power of test signal (Eq. S3)
rec_data1.run(np.max, "t")
rec_data1 *= (lambda_G / (2 * np.sqrt(np.log(2)))) * np.sqrt(2 * np.pi)
# }}}
# }}}
# {{{ Calculate receiver response as function of frequencies
# Make axis of finely spaced frequencies to feed to spline
rec_data1["nu_test"] -= carrier
Dnu = np.linspace(
    (carrier) - (rec_data1.getaxis("nu_test")[-1] / 2),
    (carrier) + (rec_data1.getaxis("nu_test")[-1] / 2),
    len(rec_data1.getaxis("nu_test")),
)
P_in = Pin_spline(Dnu)  # Generate spline of input powers
dig_filter = rec_data1 / P_in
# }}}
# {{{ Fit receiver response
A, omega, delta_nu, nu_test = symbols("A omega delta_nu nu_test", real=True)
func_form = A * sp.sinc((2 * np.pi * (nu_test - omega)) / (delta_nu)) ** 2
f = psd.lmfitdata(dig_filter)
f.functional_form = func_form
f.set_guess(
    A=dict(value=1.6e17, min=1e17, max=9e19),
    omega=dict(value=1, min=-100, max=100),
    delta_nu=dict(value=9e3, min=3, max=1e4),
)
f.fit()
fit = f.eval()
# }}}
with psd.figlist_var() as fl:
    fl.next("Receiver Response")
    fl.plot(for_plot(dig_filter), "o")
    fl.plot(for_plot(fit), color="red", alpha=0.5)
