"""
Spin Trapping Example
=====================

In this example, we perform quantitative ESR by fitting a sample with spin
trapping.

The data and fit here were used in:

Fu, X.; Xu, B.; Liyanage, H.; Zhang, C.; Kincaid, W. F.; Ford, A. L.;
Westbrook, L. G.; Brown, S. D.; DeMarco, T.; Hougland, J. L.; Franck, J. M.;
Hu, X. Ultrasound-Triggered Prodrug Activation via Sonochemically Induced
Cleavage of a 3,5-Dihydroxybenzyl Carbamate Scaffold. Chem. Sci. 2025,
`doi:10.1039.D5SC05710H <https://doi.org/10.1039/D5SC05710H>`__.

We account for slight **overmodulation** to get a close quantitative
fit.
A temporary independent-line model lets
:func:`pyspecProcScripts.peel_peaks` measure the four DMPO lines
directly.
Those measurements determine the initial amplitude, linewidth, center
field, and hyperfine splitting of the physically constrained
spin-trapping model, replacing hard-coded guesses.
Only this initial guess and the final fit are shown.
"""

import re
import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd
import pyspecProcScripts as psdpr
import sympy as sp
from numpy import pi
from pyspecdata import Q_
from scipy.special import jv


# {{{ this function describes overmodulation
def D(n, d):
    u = d.getaxis("B") * 2 * pi  # needs to be G/rad for units to work
    h_m = Q_(*d.get_prop("ModAmp")).to("G").magnitude
    arg = abs(h_m * u / 2)
    # equation below gives h_m/4, but I have already divided by h_m during
    # scaling, so this should be divided by h_m as well I also do not
    # divide by 4 here, because that would mess up my scaling
    retval = 1 / n * 1j ** (n - 1) * (jv(n - 1, arg) + jv(n + 1, arg))
    if (n - 1) % 2:
        retval["B":(None, 0)] *= -1
    return retval


# }}}
# {{{ load and baseline-correct the data
filename, C = ("250217_DMPO_5min_PBS", 101e-6)
d = psd.find_file(re.escape(filename), exp_type="francklab_esr/Warren")
d.set_prop("calibration_name", "220720new")
d = psdpr.QESR_apply_scalefactor(d)
d.set_prop("concentration", C)
if "harmonic" in d.dimlabels:
    d.chunk_auto("harmonic")
    d = d["harmonic", 0]["phase", 0]
# the following line is a baseline subtraction
d -= d["$B_0$":(3405, 3430)].mean()
d.set_plot_color_next()
plt.figure(filename, figsize=(8, 4))
# rename spectrum to "B" so it matches sympy, and reliquish control of color
d.rename("$B_0$", "B").set_plot_color(None)
C = d.get_prop("concentration")
# }}}

# {{{ configure transformations shared by the guess and final model
d.ift("B", shift=True)
d = psd.lmfitdata(d)


# The data transform converts the stored experimental data back to the
# field domain used for comparison (it is applied only to the u-domain
# data).
# In contrast, the model needs to apply D to simulate modulation, which
# is included in the residual transform.
@d.define_data_transform
def my_data_transform(d_local):
    d_local["B":0] *= 0.5
    d_local.ft("B")
    return d_local.real


@d.define_residual_transform
def my_residual_transform(d_local):
    d_local *= D(1, d_local)
    d_local["B":0] *= 0.5
    d_local.ft("B")
    return d_local.real


# }}}
# {{{ peel four independent lines to obtain model-calibrated
#      measurements
npeaks = 4
u = sp.symbols("B", real=True)
A_line, Bcenter_line, FWHM_line, L_vs_G_frac_line = (
    sp.symbols(f"{name}0:{npeaks}", real=True)
    for name in ("A_line", "Bcenter_line", "FWHM_line", "L_vs_G_frac_line")
)
# The Olivero-Longbothum approximation is
# FWHM ≈ 0.5346 lambda_L + sqrt((1-0.5346)^2 lambda_L^2 + f_G^2).
# Here f_G is 2*sqrt(log(2))*sigma, so the expression below solves this
# approximation for sigma at the requested total FWHM.
voigt_coeff = 0.5346
lambda_L_line = [FWHM_line[j] * L_vs_G_frac_line[j] for j in range(npeaks)]
sigma_line = [
    sp.sqrt(
        (FWHM_line[j] - voigt_coeff * lambda_L_line[j]) ** 2
        - (1 - voigt_coeff) ** 2 * lambda_L_line[j] ** 2
    )
    / (2 * sp.sqrt(sp.log(2)))
    for j in range(npeaks)
]
d.functional_form = sum(
    [
        A_line[j]
        * (-1j * 2 * pi * u)
        * sp.exp(
            1j * 2 * pi * Bcenter_line[j] * u
            - lambda_L_line[j] * sp.pi * abs(u)
            - sp.pi**2 * abs(u) ** 2 * sigma_line[j] ** 2
        )
        for j in range(npeaks)
    ]
)
# peel_peaks calibrates only the named amplitude, FWHM, and center parameters.
# Set the remaining shape-balance parameters first because their values affect
# the model-derived calibration.
d.set_guess(
    {
        str(symbol): {"value": 0.5, "min": 0, "max": 1}
        for symbol in L_vs_G_frac_line
    }
)
# peel_peaks depends on the functional form set above.
d = psdpr.peel_peaks(
    d,
    amplitude_parameters=[str(symbol) for symbol in A_line],
    linewidth_parameters=[str(symbol) for symbol in FWHM_line],
    center_parameters=[str(symbol) for symbol in Bcenter_line],
    close_threshold=10,
)
line_amplitudes = np.array(
    [d.guess_parameters[str(symbol)].value for symbol in A_line]
)
line_widths = np.array(
    [d.guess_parameters[str(symbol)].value for symbol in FWHM_line]
)
line_centers = np.array(
    [d.guess_parameters[str(symbol)].value for symbol in Bcenter_line]
)
# }}}
# {{{ translate peeled lines into the shared DMPO hyperfine model
A, FWHM, L_vs_G_frac, Bcenter, a_N, a_H = sp.symbols(
    "A FWHM L_vs_G_frac Bcenter a_N a_H", real=True
)
lambda_L = FWHM * L_vs_G_frac
sigma = sp.sqrt(
    (FWHM - voigt_coeff * lambda_L) ** 2 - (1 - voigt_coeff) ** 2 * lambda_L**2
) / (2 * sp.sqrt(sp.log(2)))
d.functional_form = sum(
    [
        A
        * (-1j * 2 * pi * u)
        * sp.exp(
            1j * 2 * pi * (Bcenter + I_N * a_N + I_H * a_H) * u
            - lambda_L * sp.pi * abs(u)
            - sp.pi**2 * abs(u) ** 2 * sigma**2
        )
        for I_N in (-1, 0, 1)
        for I_H in (-0.5, 0.5)
    ]
)
multiplicities = np.array([1, 2, 2, 1])
amplitude_guess = np.mean(line_amplitudes / multiplicities)
linewidth_guess = np.mean(line_widths)
center_guess = np.mean(line_centers)
splitting_guess = np.mean(np.diff(line_centers))
d.set_guess(
    {
        str(A): {
            "value": amplitude_guess,
            "min": 0,
            "max": 10 * amplitude_guess,
        },
        str(FWHM): {
            "value": linewidth_guess,
            "min": 0.1 * linewidth_guess,
            "max": 10 * linewidth_guess,
        },
        str(L_vs_G_frac): {"value": 0.5, "min": 0, "max": 1},
        str(Bcenter): {
            "value": center_guess,
            "min": center_guess - linewidth_guess,
            "max": center_guess + linewidth_guess,
        },
        str(a_H): {
            "value": splitting_guess,
            "min": 0.8 * splitting_guess,
            "max": 1.2 * splitting_guess,
        },
        str(a_N): {
            "value": splitting_guess,
            "min": 0.8 * splitting_guess,
            "max": 1.2 * splitting_guess,
        },
    }
)
# }}}
# {{{ plot only the experimental data, initial guess, and final fit
zoom_tuple = (3440, 3520)
experimental = d.data_transform(d.C)
initial_guess = d.set_to_guess().eval()
_, ax = plt.subplots(figsize=(8, 4))
psd.plot(
    experimental["B":zoom_tuple],
    "k",
    label="experimental data",
    alpha=0.7,
    ax=ax,
)
psd.plot(
    initial_guess["B":zoom_tuple],
    ":",
    label="peak-peeling initial guess",
    alpha=0.8,
    ax=ax,
)
d.fit(use_jacobian=False)
final_fit = d.eval()
psd.plot(
    final_fit["B":zoom_tuple],
    "r",
    label="final fit",
    alpha=0.7,
    ax=ax,
)
# }}}
print(f"parameters for {filename} {Q_(C, 'M').to('μM')}")
print(d.output())
print(
    "integrated amplitude:",
    final_fit.C.integrate("B", cumulative=True).integrate("B"),
)
ax.legend(loc="upper left", bbox_to_anchor=(1.1, 1.0))
ax.set_title(None)
ax.set_ylabel(None)
ax.set_xlabel("$B_0$ / G")
ax.figure.tight_layout()
plt.savefig("simulation.png", dpi=300)
plt.show()
