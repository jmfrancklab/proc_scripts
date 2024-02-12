from numpy import r_, pi
from matplotlib.pyplot import gca, cycler
import pyspecdata as psp
from scipy.optimize import leastsq, minimize, basinhopping
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
from scipy.special import erf
import sympy as sp
import numpy as np
import logging


def fit_envelope(
    s,
    min_l=10,
    threshold=0.10,
    full_width=8.434,
    direct="t2",
    show_expanding_envelope=False,
    fl=None,
):
    """
    Fit the envelope of the signal using a folded normal
    distribution to determine the lorentzian L that is 
    later fed to the L2G function for apodization
    Parameters
    ==========
    min_l:  int
            minimal possible lambda
    threshold:  float
                multiplier used in deciding the interpolation
                bounds
    full_width: float
                part of the aspect ratio for one of the diagnostic
                plots
    direct:     str
                direct dimension
    show_expanding_envelope:    boolean
                                whether the diagnostics are shown
    Returns
    =======
    env_out['lambda_l']:    float
                            Full width half max value of the fitted
                            echo envelope
    """        
    assert not s.get_ft_prop(direct), "s *must* be in time domian"
    envelope = abs(s[direct:(0, None)]).mean_all_but(direct)
    envelope = psp.lmfitdata(envelope)
    # {{{ copy/paste code for envelope
    A, lL, sigma, t = sp.symbols("A lambda_L sigma t2")
    y = A * sp.exp(-sp.pi * lL * abs(t))
    envelope.functional_form = sigma * sp.sqrt(2 / sp.pi) * sp.exp(
        -(y**2) / 2 / sigma**2
    ) + y * sp.erf(y / sp.sqrt(2 * sigma**2))
    envelope.set_guess(
        A=envelope.data.max(),
        sigma=envelope[direct, -100:].data.mean(),
        lambda_L=1 / 10e-3 / pi,
    )
    if fl:
        envelope.settoguess()
        orig_guess = envelope.eval()
    envelope.fit()
    new_guess = envelope.output()
    if fl:
        fl.push_marker()
        fl.next("envelope")
        gca().set_prop_cycle(
            cycler(alpha=[0.1, 1] + [0.5] * 5)
            + cycler(color=["g", "k", "k", "k", "r", "g", "b"])
            + cycler(ls=["-", "-", "--", ":", "-", "-", "-"])
        )
        fl.plot(orig_guess / new_guess["A"], label="guess")
        fl.plot(
            envelope / new_guess["A"],
            lw=1,
            label="signal envelope",
        )
        fl.plot(
            envelope.eval() / new_guess["A"], lw=1.1, label="least squares fit"
        )
    lsq_lambda = new_guess["lambda_L"]
    lw_range = r_[min_l : new_guess["lambda_L"] : 50j]
    amount_over = np.zeros_like(lw_range)
    # where does A exp(-π λ t) decay to 2σ?
    # at -ln(A/2σ)/π λ
    # = ln(2σ/A) / π λ
    t_at_exp_end = (
        np.log(new_guess["A"] / 2 / new_guess["sigma"]) / pi / lsq_lambda
    )
    for j, newL in enumerate(lw_range):
        new_guess.update(lambda_L=newL)
        envelope.set_guess(new_guess)
        envelope.settoguess()
        points_over = (
            envelope[direct:(0, t_at_exp_end)]
            - envelope.eval()[direct:(0, t_at_exp_end)]
        )
        points_over[lambda x: x < 0] = 0
        points_over.run(lambda x: np.sqrt(abs(x) ** 2)).mean()
        amount_over[j] = points_over.item()
    # }}}
    l = "$\\lambda_L$"
    env_expansion = psp.nddata(amount_over / amount_over.max(), [-1], [l])
    env_expansion.setaxis(l, lw_range).set_units(l, "Hz")
    env_expansion.name("norm of points\noutside envelope")
    if fl and show_expanding_envelope:
        fl.next("expanding envelope", figsize=r_[1, 0.3] * full_width)
        fl.plot(env_expansion)
    norm_max = env_expansion.max().item().real
    norm_min = env_expansion.min().item().real
    opt_lambda = env_expansion.invinterp(
        l, norm_min * (1 - threshold) + norm_max * threshold, kind="linear"
    )
    #logging.debug(
    #    "opt_lambda",
    #    opt_lambda,
    #    "at",
    #    norm_min * (1 - threshold) + norm_max * threshold,
    #    "out of",
    #    lw_range,
    #)
    if fl and show_expanding_envelope:
        fl.plot(opt_lambda, "o")
    new_guess.update(lambda_L=opt_lambda.getaxis(l).item().real)
    envelope.set_guess(new_guess)
    envelope.settoguess()
    if fl:
        fl.next("envelope")
    env_out = envelope.output()
    if fl:
        fl.plot(
            envelope.eval() / env_out["A"], lw=1.1, label=r"optimal envelope"
        )
    t2 = envelope.fromaxis(direct)
    # lsq
    lsq = np.exp(-pi * lsq_lambda * abs(t2))
    if fl:
        fl.plot(
            L2G(env_out["lambda_L"], criterion="energy")(s.fromaxis(direct)),
            label = 'apodization function')
        fl.pop_marker()
    return env_out["lambda_L"]


def L2G(
    lambda_L,
    criterion="energy",
):
    """
    Parameters
    ==========
    lambda_L:   float
                Full width half max of the signal envelope
    criterion:  str
                the final function will be made in order to 
                generate signal either with equal energy or
                with equal linewidth
    Returns
    =======
    Apodization function
    """

    assert np.isscalar(lambda_L)
    if criterion == "energy":
        # equal energy:
        return lambda t2: np.exp(
            -0.5 * pi**3 * lambda_L**2 * t2**2 + pi * lambda_L * abs(t2)
        )
    elif criterion == "width":
        # equal linewidth
        return lambda t2: np.exp(
            pi
            * lambda_L
            * (-pi * lambda_L * t2**2 / 4 / np.log(2) + abs(t2))
        )
