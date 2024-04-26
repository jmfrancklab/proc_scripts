from pint import UnitRegistry
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import QESR_scalefactor
from sympy import symbols, Symbol, latex
from sympy import exp as s_exp
import numpy as np

logger = init_logging(level="debug")
fieldaxis = "$B_0$"
exp_type = "francklab_esr/alex"
with figlist_var() as fl:
    for filenum, (thisfile, calibration, diameter, background) in enumerate(
        [
            (
                "230418_200uM_TEMPOL_3_redo2.DSC",
                "230202",
                "QESR caps",
                find_file("230418_water.DSC", exp_type=exp_type)["harmonic", 0],
            ),
        ]
    ):
        # {{{ load data in and rescale and center about 0
        d = find_file(thisfile, exp_type=exp_type)
        if "harmonic" in d.dimlabels:
            d = d["harmonic", 0]
        color = d.get_plot_color()
        d /= QESR_scalefactor(d, calibration_name=calibration, diameter_name=diameter)
        #center_field = d.getaxis(fieldaxis)[r_[0, -1]].mean()
        d -= d[fieldaxis, -100:].data.mean()
        #d.setaxis(fieldaxis, lambda x: x - center_field)
        background = ndshape(d).alloc()
        background.copy_axes(d)
        d -= background
        # }}}
        # {{{ Plot starting spectra
        d.set_units(fieldaxis, "G")
        fl.next("B domain")
        fl.plot(d, label="actual data")
        d_axis = (d[fieldaxis][:])
        d.ift(fieldaxis,shift = True)
        fl.next("u domain")
        fl.plot(d, label="real data")
        fl.plot(d.imag, label="imag data")
        # }}}
        # {{{ Make list of guesses
        centers = [
            3.4894e3,
            3.5064e3,
            3.5235e3,
        ]  # make a list of the centers of the starting data for the fits
        lambda_l_guess = [1.2, 1.2, 1.2]
        lambda_g_guess = [1.2, 1.2, 1.2]
        amp_guess = [80, 80, 80]
        # }}}
        # {{{Initiate lmfit
        d.rename(fieldaxis, "u")
        u = d.fromaxis("u")
        f = lmfitdata(d)
        # }}}
        # {{{ make symbols and prep lists for symbols
        all_symb = {symbname:[symbols("%s%d"%(symbname,(j+1)),real=True) for j in range(3)] for symbname in ["lambda_l","lambda_g","amp","B"]}
        u = symbols("u", real=True)
        # }}}
        # {{{ index and make the functions using the symbols - we will have 3 lorentzians, 3 gaussians and 3 shifts
        Lorentz = []  # the list needs place holders for my loop to work
        Gauss = []
        Shifts = []
        for j in range(3):
            Lorentz.append(s_exp(-pi * all_symb["lambda_l"][j] * abs(u)))
            Gauss.append(
                s_exp(-(pi**2 * all_symb["lambda_g"][j] ** 2 * u**2) / (4 * np.log(2)))
            )
            Shifts.append(s_exp(+1j * pi * 2 * u * all_symb["B"][j]))
        # }}}
        # {{{ make a fitting function as a sum of the 3 voigt lines using the function lists we just made
        deriv = -1j * 2 * pi * u
        myfit = 0
        for j in range(len(Gauss)):
            myfit += all_symb["amp"][j] * Gauss[j] * Lorentz[j] * Shifts[j] * deriv
        # }}}
        f.functional_form = myfit
        # {{{ set your guesses
        all_lambda_g = {
            "lambda_g%d" % (j + 1): {"value": lambda_g_guess[j], "min": 0.1, "max": 3}
            for j in range(3)
        }
        all_amp = {
            "amp%d" % (j + 1): {"value": amp_guess[j], "min": 5e-6, "max": 5e6}
            for j in range(3)
        }
        all_lambda_l = {
            "lambda_l%d" % (j + 1): {"value": lambda_l_guess[j], "min": 0.1, "max": 3}
            for j in range(3)
        }
        all_shift = {
            "B%d" % (j + 1): {"value": centers[j], "min": 3.4604e3, "max": 3.5800e3}
            for j in range(3)
        }
        # }}}
        # {{{plot your guess
        f.set_guess(
            **all_lambda_l,
            **all_lambda_g,
            **all_shift,
            **all_amp,
        )
        f.settoguess()
        guess = f.eval()
        fl.plot(guess,label='guess')
        guess.ft("u", shift=True)
        fl.next("B domain")
        guess.set_units('u','kG')
        fl.plot(guess, label="guess")
        fl.show();quit()
        # }}}
        # {{{fit
        f.fit()
        fit = f.eval()
        # }}}
        fit.ft("u", shift=True)
        fit.ift("u")
        fl.next("u domain")
        fl.plot(fit, label="fit")
        fl.next("B domain")
        fit.ft("u")
        fit.set_units("u", "G")
        fl.plot(fit, label="fit")
