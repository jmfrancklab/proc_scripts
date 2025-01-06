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
lambda_G = 0.4e3  # Width for Gaussian convolution
# {{{ Calculate input power (acquired on Oscilloscope)
#     To determine the input power, just take a capture of the signal in the
#     time domain and fit to a complex exponential
# There are less nodes acquired for the control case so I need to separately
# define them
file1_nodes = psd.find_file(
    re.escape(file1),
    exp_type="ODNP_NMR_comp/noise_tests",
    return_list=True,
)
file2_nodes = psd.find_file(
    re.escape(file2),
    exp_type="ODNP_NMR_comp/noise_tests",
    return_list=True,
)
# {{{ Sort node names based on frequency output
#     by source
frqs_str_kHz = sorted(file1_nodes, key=lambda x: int(x.split("_")[1]))
# }}}
# {{{ Make empty nddata to drop the calculated $V^{2}$ into with corresponding
#     frequency ($\nu$) output by AFG source, and set the frequencies based on
#     the node names
control = psd.ndshape([len(frqs_str_kHz)], ["nu"]).alloc()
control.setaxis(
    "nu",
    np.array(list(int(j.split("_")[1]) * 1e3 for j in frqs_str_kHz)),
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
    control["nu", j] = V_amp**2
# {{{ make spline for power going into RX box
Pin_spline = control.spline_lambda()
# }}}
# }}}
# {{{ Calculate $dg^2$
frq_kHz = np.array([int(j.split("_")[1]) for j in file2_nodes])
rec_data = (
        psd.ndshape([len(file2_nodes)],
        ["nu_test"]).alloc().set_units("nu_test","Hz")
        )
rec_data["nu_test"] = frq_kHz
# Lists for calculated dg^2, frequency output by source, and
# the offset of the output frequency and $\nu_{RX,LO}$
for j, nodename in enumerate(file2_nodes):
    d = psd.find_file(
        file2,
        exp_type="ODNP_NMR_comp/noise_tests",
        expno=nodename,
    )
    carrier = (
        d.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
    )  # $\nu_{RX,LO}$
    SW = str(d.get_prop("acq_params")["SW_kHz"])
    d.set_units("t", "s")
    d.rename(
        "nScans", "capture"
    )  # To be more consistent with the oscilloscope data rename
    #    the nScans dimension
    acq_time = np.diff(d.getaxis("t")[np.r_[0, -1]]).item()
    d.ft("t", shift=True)  # $dg\sqrt{s/Hz}$
    d = abs(d) ** 2  # $dg^{2}*s/Hz$
    d.mean("capture")
    d /= acq_time  # $dg^{2}/Hz$
    # Convolve using $\lambda_{G}$ specified above
    d.convolve("t", lambda_G, enforce_causality=False)
    d.ift("t")
    d.run(np.conj)  # Empirically needed to give offset
    #                        that increases with field
    d.ft("t")
    d["t"] += carrier
    # }}}
    # {{{ Calculate power of test signal (Eq. S3)
    d.run(np.max, "t")
    d *= (lambda_G / (2 * np.sqrt(np.log(2)))) * np.sqrt(np.pi)
    # }}}
    rec_data["nu_test",j] = d
# set x axis to $\Delta\nu$ and rename accordingly
rec_data.rename("nu_test", "nu_offset")
print(carrier)
print(rec_data["nu_offset"])
rec_data["nu_offset"] -= int(carrier)
Dnu = np.linspace(
    (carrier) - (rec_data.getaxis("nu_offset")[-1] / 2),
    (carrier) + (rec_data.getaxis("nu_offset")[-1] / 2),
    len(rec_data.getaxis("nu_offset")),
)
print(Dnu)
# }}}
# {{{ Calculate receiver response as function of frequencies
P_in = Pin_spline(Dnu)  # Generate spline of input powers
dig_filter = rec_data / P_in  # $dg^2/V_{in}^2$
# }}}
with psd.figlist_var() as fl:
    fl.next("P in")
    P_in.human_units(scale_data=True)
    fl.plot(P_in)
    fl.next("test")
    dig_filter.human_units(scale_data=True)
    fl.plot(dig_filter,'o')
    fl.show();quit()
# {{{ Fit receiver response
A, omega, delta_nu, nu_offset = symbols(
    "A omega delta_nu nu_offset", real=True
)
func_form = A * sp.sinc((2 * np.pi * (nu_offset - omega)) / (delta_nu)) ** 2
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
    plt.ylabel(r"$\mathrm{dg_{%s\ \mathrm{kHz}}}/ \mu V$" % SW)
