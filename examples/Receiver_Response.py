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
from numpy import r_, pi
import pyspecdata as psd
from sympy import symbols
import sympy as sp
import re

Dnu_name = r"$\Delta\nu$"
lambda_G = 0.4e3  # Width for Gaussian convolution
data_dir = "ODNP_NMR_comp/noise_tests"
file1 = "240123_10mV_AFG_GDS_5mV_100MSPS_analytic.h5"
file2 = "240117_afg_sc_10mV_3p9kHz_zoom.h5"
# There are less nodes acquired for the control case (since we assume it's
# relatively flat) so I need to separately define them
file1_nodes = psd.find_file(
    re.escape(file1),
    exp_type=data_dir,
    return_list=True,
)
control_frqs = np.array([int(j.split("_")[1]) for j in file1_nodes])*1e3
file2_nodes = psd.find_file(
    re.escape(file2),
    exp_type=data_dir,
    return_list=True,
)
# $\nu_{test}$
nu_test = np.array([int(j.split("_")[1]) for j in file2_nodes])
# {{{ Calculate input V (acquired on Oscilloscope)
# {{{ Make empty nddata to drop the calculated V into with corresponding
#     frequency ($\nu$) output by AFG source, and set the frequencies based on
#     the node names
control = psd.ndshape([len(control_frqs)], ["nu_test"]).alloc()
control.setaxis("nu_test", control_frqs).set_units("nu_test", "Hz")
# }}}
for j, nodename in enumerate(file1_nodes):
    rf_frq = control["nu_test"][j]
    d = psd.find_file(
        file1,
        expno=nodename,
        exp_type=data_dir,
    )
    # {{{ fit signal in t domain to complex exponential
    A, omega, phi, t = symbols("A omega phi t", real=True)
    f = psd.lmfitdata(d)
    f.functional_form = A * sp.exp(1j * 2 * pi * omega * t + 1j * phi)
    f.set_guess(
        A=dict(value=5e-3, min=1e-4, max=1),
        omega=dict(
            value=rf_frq,
            min=rf_frq - 1e4,
            max=rf_frq + 1e4,
        ),
        phi=dict(value=0.75, min=-pi, max=pi),
    )
    f.fit(use_jacobian=False)
    f.eval()
    # }}}
    V_amp = f.output("A")
    control["nu_test", j] = abs(V_amp) * 1e6 # μV
# {{{ make spline for power going into RX box
control.sort("nu_test")
control.rename("nu_test",Dnu_name) # since we will be applying $\Delta\nu$
#                                    axis to spline
Pin_spline = control.spline_lambda()
# }}}
# }}}
# {{{ Calculate dg
for j, nodename in enumerate(file2_nodes):
    d = psd.find_file(
        file2,
        exp_type=data_dir,
        expno=nodename,
    )
    carrier = (
        d.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
    )  # $\nu_{RX,LO}$
    if j == 0:
        # Allocate and array that's shaped like one of the datasets
        # acquired (with 100 captures and a direct t dimension)
        # but add axis to store the frequency of the test signal
        rec_data = (d.shape + ("nu_test", len(nu_test))).alloc()
        rec_data.setaxis("t", d.getaxis("t")).set_units("t", "s")
        rec_data.setaxis("nu_test", nu_test).set_units("nu_test", "Hz")
    rec_data["nu_test", j] = d
    SW = str(d.get_prop("acq_params")["SW_kHz"]) # For labeling final plot
rec_data.sort("nu_test")
rec_data.rename("nScans", "capture")  # To be more consistent with the
#                                        oscilloscope data rename the nScans
#                                        dimension
acq_time = np.diff(rec_data["t"][r_[0, -1]]).item()
# {{{ Calculate PSD for each frequency (we will calculate power from the A
#     of the convolved test signal)
rec_data.run(np.conj)  # Empirically needed to give offset
#                        that increases with field
rec_data.ft("t", shift=True)  # $dg\sqrt{s/Hz}$
rec_data = abs(rec_data) ** 2  # $dg^{2}*s/Hz$
rec_data.mean("capture")
rec_data /= acq_time  # $dg^{2}/Hz$
# Convolve using $\lambda_{G}$ specified above
rec_data.convolve("t", lambda_G, enforce_causality=False)
# }}}
# {{{ Calculate power of test signal (Eq. S3)
rec_data["t"] += carrier
rec_data.run(np.max, "t") # Takes maximum of PSD $\rightarrow dg^{2}$
rec_data *= (lambda_G / (2 * np.sqrt(np.log(2)))) * np.sqrt(pi)
rec_data.run(np.sqrt) # dg
# }}}
# }}}
# {{{ Calculate receiver response as function of frequencies
# Make axis of finely spaced frequencies to feed to spline
rec_data["nu_test"] = rec_data["nu_test"] - carrier # center data at 0 MHz
#                                                     thus converting to $\Delta\nu$
#                                                     rather than $\nu_{test}$
rec_data.rename("nu_test",Dnu_name)
Dnu = np.linspace(
    (carrier) - (rec_data.getaxis(Dnu_name)[-1] / 2),
    (carrier) + (rec_data.getaxis(Dnu_name)[-1] / 2),
    len(rec_data.getaxis(Dnu_name)),
)
control = Pin_spline(Dnu)  # Generate spline of input powers
dig_filter = rec_data / control # dg/μV
dig_filter.name(r"$\mathrm{dg_{%s\ \mathrm{kHz}}}/ \mu V$" % SW)
# }}}
# {{{ Fit receiver response
A, omega, delta_nu, Dnu_name = symbols(r"A omega delta_nu $\Delta\nu$", real=True)
f = psd.lmfitdata(dig_filter)
f.functional_form = A * abs(sp.sinc((2 * pi * (Dnu_name - omega)) / (delta_nu)))
f.set_guess(
    A=dict(value=400, min=1, max=9e3),
    omega=dict(value=1, min=-100, max=100),
    delta_nu=dict(value=9e3, min=3, max=1e4),
)
f.fit()
fit = f.eval()
# }}}
with psd.figlist_var() as fl:
    fl.next("Receiver Response")
    fl.plot(dig_filter, "o")
    fl.plot(fit, color="red", alpha=0.5)
