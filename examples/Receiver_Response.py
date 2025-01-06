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
import pylab as plt
import pyspecdata as psd
from sympy import symbols
from scipy.interpolate import CubicSpline
import sympy as sp
import re

file1 = "240123_10mV_AFG_GDS_5mV_100MSPS_analytic.h5"
file2 = "240117_afg_sc_10mV_3p9kHz_zoom.h5"
# {{{ Calculate input power (acquired on Oscilloscope)
#     To determine the input power, just take a capture of
#     the signal in the time domain and fit to a complex
#     exponential
all_nodes = psd.find_file(
    re.escape(file1),
    exp_type="ODNP_NMR_comp/noise_tests",
    return_list=True,
)
# {{{ Sort node names based on frequency output
#     by source
frqs_str_kHz = sorted(
    all_nodes, key=lambda x: int(x.split("_")[1])
)
# }}}
# {{{ Make empty nddata to drop the calculated
#     $V^{2}$ into with corresponding frequency
#     ($\nu$) output by AFG source, and set the
#     frequencies based on the node names
control = psd.ndshape([len(frqs_str_kHz)], ["nu"]).alloc()
control.setaxis(
    "nu",
    np.array(
        list(
            int(j.split("_")[1]) * 1e3 for j in frqs_str_kHz
        )
    ),
).set_units("nu", "Hz")
# }}}
for j, nodename in enumerate(frqs_str_kHz):
    rf_frq = control["nu"][j]
    d = psd.find_file(
        re.escape(file1),
        expno=nodename,
        exp_type="ODNP_NMR_comp/noise_tests",
    )
    # {{{ fit signal in t domain to complex exponential
    A, omega, phi, t = symbols("A omega phi t", real=True)
    f = psd.lmfitdata(d)
    f.functional_form = A * sp.exp(
        1j * 2 * np.pi * omega * t + 1j * phi
    )
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
    control["nu", j] = V_amp**2
# {{{ make spline for power going into RX box
t_Pin_cs = CubicSpline(
    np.array(control.getaxis("nu")), control.data
)
# }}}
# }}}
# {{{ Calculate $dg^2$
lambda_G = 0.4e3  # Width for Gaussian convolution
all_nodes = psd.find_file(
    re.escape(file2),
    exp_type="ODNP_NMR_comp/noise_tests",
    return_list=True,
)
# Lists for calculated dg^2, frequency output by source, and
# the offset of the output frequency and $\nu_{RX,LO}$
rec_data = []
nu_test_kHz = []  # $\nu_{test}$
for j, nodename in enumerate(all_nodes):
    d = psd.find_file(
        file2,
        exp_type="ODNP_NMR_comp/noise_tests",
        expno=nodename,
    )
    carrier = (
        d.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
    )  # $\nu_{RX,LO}$
    rec_data.append(d)
    SW = str(d.get_prop("acq_params")["SW_kHz"])
    nu_test_kHz.append(d.get_prop("afg_frq"))
nu_test_kHz.sort()
# Make additional dimension to store the frequencies of test
# signal
rec_data = psd.concat(rec_data, "nu_test")
rec_data["nu_test"] = np.array(nu_test_kHz) * 1e3
# }}}
rec_data.sort("nu_test")
rec_data.set_units("t", "s")
rec_data.set_units("nu_test", "s")
rec_data.rename(
    "nScans", "capture"
)  # To be more consistent with the oscilloscope data rename
#    the nScans dimension
acq_time = np.diff(
    rec_data.getaxis("t")[np.r_[0, -1]]
).item()
rec_data.ft("t", shift=True)  # $dg\sqrt{s/Hz}$
rec_data = abs(rec_data) ** 2  # $dg^{2}*s/Hz$
rec_data.mean("capture")
rec_data /= acq_time  # $dg^{2}/Hz$
# Convolve using $\lambda_{G}$ specified above
rec_data.convolve("t", lambda_G, enforce_causality=False)
rec_data.ift("t")
rec_data.run(np.conj)  # Empirically needed to give offset
#                        that increases with field
rec_data.ft("t")
rec_data["t"] += carrier
# }}}
# {{{ Calculate power of test signal (Eq. S3)
rec_data.run(np.max, "t")
rec_data *= (lambda_G / (2 * np.sqrt(np.log(2)))) * np.sqrt(
    np.pi
)
# }}}
# set x axis to $\Delta\nu$ and rename accordingly
rec_data.rename("nu_test", "nu_offset")
rec_data["nu_offset"] -= carrier
Dnu = np.linspace(
    (carrier / 1e3)
    - (rec_data.getaxis("nu_offset")[-1] / 2),
    (carrier / 1e3)
    + (rec_data.getaxis("nu_offset")[-1] / 2),
    len(rec_data.getaxis("nu_offset")),
)
# }}}
# {{{ Calculate receiver response as function of frequencies
P_in = t_Pin_cs(Dnu)  # Generate spline of input powers
dig_filter = rec_data / P_in  # $dg^2/V_{in}^2$
print(dig_filter["nu_offset"])
# }}}
# {{{ Fit receiver response
A, omega, delta_nu, nu_offset = symbols(
    "A omega delta_nu nu_offset", real=True
)
func_form = (
    A
    * sp.sinc(
        (2 * np.pi * (nu_offset - omega)) / (delta_nu)
    )
    ** 2
)
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
# Take sqrt for units of $dg/V$
dig_filter.data = np.sqrt(dig_filter.data)
fit.data = np.sqrt(fit.data)
dig_filter /= 1e6  # Units of $dg/\mu V$
fit /= 1e6  # Units of $dg/\mu V$
with psd.figlist_var() as fl:
    fl.next("Receiver Response")
    fl.plot(dig_filter, "o")
    fl.plot(fit, color="red", alpha=0.5)
    plt.xlabel(r"$\Delta \nu$ / kHz")
    plt.ylabel(
        r"$\mathrm{dg_{%s\ \mathrm{kHz}}}/ \mu V$" % SW
    )
