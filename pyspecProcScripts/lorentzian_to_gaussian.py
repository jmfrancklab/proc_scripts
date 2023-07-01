from numpy import r_, pi
import pyspecdata as psp
from pyspecdata import nddata
from scipy.optimize import leastsq, minimize, basinhopping
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
from scipy.special import erf
import sympy as sp
from numpy import amin


def L2G(
    s, criterion="energy", min_l=10, threshold=0.10, full_width=8.434, fl=None
):
    assert not s.get_ft_prop("t2")
    # s *must* be in time domian
    envelope = abs(s["t2":(0, None)]).mean_all_but("t2")
    envelope = lmfitdata(envelope)
    # {{{ copy/paste code for envelope
    A, lL, sigma, t = sp.symbols("A lambda_L sigma t2")
    y = A * sp.exp(-sp.pi * lL * abs(t))
    envelope.functional_form = sigma * sp.sqrt(2 / sp.pi) * sp.exp(
        -(y**2) / 2 / sigma**2
    ) + y * sp.erf(y / sp.sqrt(2 * sigma**2))
    envelope.set_guess(
        A=envelope.data.max(),
        sigma=envelope["t2", -100:].data.mean(),
        lambda_L=1 / 10e-3 / pi,
    )
    envelope.settoguess()
    show_guess = False
    if show_guess:
        fl.plot(envelope.eval(), label="guess")
    envelope.fit()
    new_guess = envelope.output()
    fl.next("envelope")
    fl.plot(
        envelope / new_guess["A"], "k", alpha=1, lw=1, label="signal envelope"
    )
    gca().set_prop_cycle(
        cycler(alpha=[0.5]) * cycler(color=["k", "k", "r", "g", "b"])
        + cycler(ls=["--", ":", "-", "-", "-"])
    )
    fl.plot(
        envelope.eval() / new_guess["A"], lw=1.1, label="least squares fit"
    )
    lsq_lambda = new_guess["lambda_L"]
    # lw_range = r_[new_guess["lambda_L"] : new_guess["lambda_L"] / 3 : 50j]
    lw_range = r_[min_l : new_guess["lambda_L"] : 50j]
    amount_over = zeros_like(lw_range)
    for j, newL in enumerate(lw_range):
        new_guess.update(lambda_L=newL)
        envelope.set_guess(new_guess)
        envelope.settoguess()
        # if j == 0:
        #    fl.plot(envelope.eval()/new_guess['A'], label=r"min $\lambda_L$")
        points_over = envelope - envelope.eval()
        points_over[lambda x: x < 0] = 0
        points_over.run(lambda x: sqrt(abs(x) ** 2)).mean()
        amount_over[j] = points_over.item()
    # }}}
    # fl.plot(s_noise/new_guess['A'], 'k', alpha=0.1, lw=1, label="noise (inactive $\\Delta p$)")
    fl.next("expanding envelope", figsize=r_[1, 0.3] * full_width)
    l = "$\\lambda_L$"
    env_expansion = nddata(amount_over / amount_over.max(), [-1], [l])
    env_expansion.setaxis(l, lw_range).set_units(l, "Hz")
    env_expansion.name("norm of points\noutside envelope")
    fl.plot(env_expansion)
    norm_max = env_expansion.max().item().real
    norm_min = env_expansion.min().item().real
    opt_lambda = env_expansion.invinterp(
        l, norm_min * (1 - threshold) + norm_max * threshold, kind="linear"
    )
    print(
        "opt_lambda",
        opt_lambda,
        "at",
        norm_min * (1 - threshold) + norm_max * threshold,
        "out of",
        lw_range,
    )
    fl.plot(opt_lambda, "o")
    new_guess.update(lambda_L=opt_lambda.getaxis(l).item().real)
    envelope.set_guess(new_guess)
    envelope.settoguess()
    fl.next("envelope")
    env_out = envelope.output()
    fl.plot(envelope.eval() / env_out["A"], lw=1.1, label=r"optimal envelope")
    t2 = envelope.fromaxis("t2")
    # lsq
    lsq = exp(-pi * lsq_lambda * abs(t2))
    # matched:
    matched = exp(-pi * env_out["lambda_L"] * abs(t2))
    if criterion == "energy":
        # equal energy:
        L2G = exp(
            -0.5 * pi**3 * env_out["lambda_L"] ** 2 * t2**2
            + pi * env_out["lambda_L"] * abs(t2)
        )
    elif criterion == "width":
        # equal linewidth
        L2G = exp(
            pi
            * env_out["lambda_L"]
            * (-pi * env_out["lambda_L"] * t2**2 / 4 / log(2) + abs(t2))
        )
    fl.plot(L2G)
    return L2G
