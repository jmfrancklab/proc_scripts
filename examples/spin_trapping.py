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

The strategy here involves two key elements:

    - We account for a slight **overmodulation** to get a very close fit to the
      spectrum.  This is also important for precise quantification.
    - We perform the fit **stepwise**, so that each step is not too far off and
      therefore guaranteed to converge.
      It's not necessary to show what the result of each step,
      but they are shown here so the reader can see how the fit is done.
"""

import matplotlib.pyplot as plt
import pyspecProcScripts as psdpr
import re
from pyspecdata import Q_
from numpy import pi
import sympy as sp
from scipy.special import jv
import pyspecdata as psd


# {{{ this function is used to describe overmodulation.  It's OK to not fully
#     understand where this comes from, but I can provide references.
def D(n, d):
    u = d.getaxis("B") * 2 * pi  # needs to be G/rad for units to work
    h_m = Q_(*d.get_prop("ModAmp")).to("G").magnitude
    mod_calib = 1
    h_m *= mod_calib
    arg = abs(h_m * u / 2)
    # equation below gives h_m/4, but I have already divided by h_m during
    # scaling, so this should be divided by h_m as well I also do not
    # divide by 4 here, because that would mess up my scaling
    retval = 1 / (n) * 1j ** (n - 1) * (jv(n - 1, arg) + jv(n + 1, arg))
    if (n - 1) % 2:
        retval["B":(None, 0)] *= -1
    return retval


# }}}
# {{{ load the data, and plot the raw data
filename, C = ("250217_DMPO_5min_PBS", 101e-6)
d = psd.find_file(re.escape(filename), exp_type="francklab_esr/Warren")
d.set_prop("calibration_name", "220720new")
d = psdpr.QESR_apply_scalefactor(d)
d.set_prop("concentration", C)
print(d.get_prop("ModAmp"))
if "harmonic" in d.dimlabels:
    d.chunk_auto("harmonic")
    d = d["harmonic", 0]["phase", 0]
# the following line is a baseline subtraction
d -= d["$B_0$":(3405, 3430)].mean()
d.set_plot_color_next()
plt.figure(filename, figsize=(8, 4))
# rename spectrum to "B" so it matches sympy, and reliquish control of color
d.rename("$B_0$", "B").set_plot_color(None)
# }}}
# {{{ the variables we need for our fit
C = d.get_prop("concentration")
A, u, l_L, Bcenter, sigma, a_N, a_H = sp.symbols(
    "A B lambda_L Bcenter sigma a_N a_H", real=True
)
d.ift("B", shift=True)  # Inverse fourier transform our spectral d,
#                         since it's easier to define the lineshapes in the
#                         "time" domain vs. the field/frequency domain
d = psd.lmfitdata(d)  # converts the "d" object to a form that
#                       contains all the functions needed for fitting
thefunction = 0
the_guesses = {}
# }}}
# {{{ Construct the function that represents the spectrum, which is a sum of
#     several Voigt lines
#
# loop over all our peaks, where I_N and I_H give the nuclear spin quantum
# numbers of the hyperfine coupled H and N (NOTE: for TEMPOL, there is no
# H, only N).
#
# Inside the loop, we define the functional form that determines the lineshape.
#
# To get a decent fit, we also need a somewhat reasonable guess for all
# the parameters (linewidth λ/l_L, amplitude A, center field Bcenter) that
# we use to in the function.
for I_N in (-1, 0, 1):
    for I_H in (-0.5, 0.5):
        thefunction += (
            A
            * (-1j * 2 * pi * u)  # here "u" is the Fourier conjugate of "B"
            #                     ("u" is to "B" as time is to frequency)
            * sp.exp(1j * 2 * pi * (Bcenter + I_N * a_N + I_H * a_H) * u)
            * sp.exp(-l_L * sp.pi * abs(u) - sp.pi**2 * abs(u) ** 2 * sigma**2)
        )
d.functional_form = thefunction
# }}}
# {{{ Guess the quantitites (hyperfine, concentration, etc) that we need to
#     describe our spectrum start out with all guesses the same
the_guesses.update(
    {
        str(A): {
            # peak amplitude
            "value": 2e4 * C,
            "min": 0,
            "max": 20e4 * C,
        },
        str(l_L): {
            # Lorentzian broadening
            "value": 0.01,
            "min": 0.005,
            "max": 1,
        },
        str(sigma): {
            # Gaussian broadening
            "value": 0.5,
            "min": 0.3,
            "max": 1,
        },
        str(Bcenter): {
            # Center field
            "value": 3482,
            "min": 3480,
            "max": 3484,
        },
        str(a_H): {
            # Hyperfine splitting for hydrogen
            # (Note you can remove this, and the associated term inside the
            # exponential for a nitroxide)
            "value": 15,
            "min": 14.75,
            "max": 16,
        },
        str(a_N): {
            "value": 15,
            "min": 14.75,
            "max": 16,
        },
    }
)
# }}}
d.set_guess(the_guesses)


# {{{ It's better to define our model equations in the u-domain (in Fourier
#     space), but we need to define functions that transform us back into the
#     Field (B)-domain.  That's what these are.
#
#     ALSO: the data includes overmodulation, so we move it back to the
#     B-domain with a standard FT, while our model equations do not include
#     overmodulation, so we need to include overmodulation when moving the
#     model equation into the B-domain
@d.define_data_transform
def my_data_transform(d):
    d["B":0] *= 0.5
    d.ft("B")
    return d.real


@d.define_residual_transform
def my_residual_transform(d):
    d *= D(1, d)
    d["B":0] *= 0.5
    d.ft("B")
    return d.real


# }}}
ax = plt.gca()
# {{{ plot the data
#
# Above, I inverse Fourier transformed the d, so to get it back in
# the field domain (which is equivalent to the frequency domain), I need
# to make a copy (.C) and Fourier Transform (ft) it here:
thedata = d.C.ft("B")
zoom_tuple = (3440, 3520)  # the region we want to plot
psd.plot(
    thedata["B":zoom_tuple], "k", label="experimental data", alpha=0.7, ax=ax
)
# }}}
# {{{ show the initial guess
forplot = d.settoguess().eval()
psd.plot(forplot["B":zoom_tuple], label="initial guess", alpha=0.7, ax=ax)
# }}}
# {{{ vary only the center field
for name, par in d.guess_parameters.items():
    if not name.startswith("B"):
        par.vary = False
print("about to fit field")
d.fit(use_jacobian=False)
forplot = d.eval()
psd.plot(forplot["B":zoom_tuple], label="vary only field", alpha=0.7, ax=ax)
# }}}
d.guess_parameters = d.fit_parameters
# {{{ add variation of center field and splitting
for name, par in d.guess_parameters.items():
    if name.startswith("a_"):
        par.vary = True
print("about to fit splitting")
d.fit(use_jacobian=False)
forplot = d.eval()
psd.plot(
    forplot["B":zoom_tuple], label="vary splitting, as well", alpha=0.7, ax=ax
)
# }}}
d.guess_parameters = d.fit_parameters
# {{{ vary splitting, center field, and amplitude
for name, par in d.guess_parameters.items():
    if name.startswith("A"):
        par.vary = True
print("about to fit amplitude")
d.fit(use_jacobian=False)
forplot = d.eval()
psd.plot(
    forplot["B":zoom_tuple], label="vary amplitude, as well", alpha=0.7, ax=ax
)
# }}}
d.guess_parameters = d.fit_parameters
# {{{ vary everything
for name, par in d.guess_parameters.items():
    par.vary = True
print("about to fit everything")
d.fit(use_jacobian=False)
thefit = d.eval()
psd.plot(thefit["B":zoom_tuple], label="final fit", alpha=0.7, ax=ax)
# }}}
print(f"parameters for {filename} {Q_(C, 'M').to('μM')}")
out = d.output()
print(out)
print(
    "integrated amplitude:",
    thefit.C.integrate("B", cumulative=True).integrate("B"),
)
# {{{ final plot layout
ax.legend(loc="upper left", bbox_to_anchor=(1.1, 1.0))
plt.tight_layout()
ax = plt.gca()
ax.set_title(None)
ax.set_ylabel(None)
ax.set_xlabel("$B_0$ / G")
ax.figure.tight_layout()
plt.savefig("simulation.png", dpi=300)
plt.show()
# }}}
