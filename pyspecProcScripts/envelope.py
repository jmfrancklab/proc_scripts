from numpy import r_, pi
from matplotlib.pyplot import gca, cycler
import pyspecdata as psp
import sympy as sp
import numpy as np
import logging


def fit_envelope(
    s,
    min_l=10,
    threshold=0.10,
    full_width=8.434,
    direct="t2",
    plot_name="envelope",
    show_expanding_envelope=False,
    mult_two=False,
    fl=None,
):
    """Fit the homogeneous Lorentzian linewidth of an FID decay.

    Parameters
    ==========
    s : nddata
        Time-domain FID-side signal.
    min_l : float
        Smallest linewidth, in Hz, considered by the expanding-envelope
        refinement.
    threshold : float
        Fraction of the expanding-envelope range used to select the refined
        linewidth.
    full_width : float
        Figure width used for the optional expanding-envelope diagnostic.
    direct : str
        Direct time dimension.
    plot_name : str
        Figure name for the envelope-fit diagnostic.
    show_expanding_envelope : bool
        Plot the expanding-envelope refinement when diagnostics are enabled.
    mult_two : bool
        Restore a half-weighted zero-time point before fitting. Use this for
        data returned by :func:`fid_side_from_echo`.
    fl : figlist or None
        Optional diagnostic figure list.

    Returns
    =======
    homogeneous_linewidth : float
        Positive Lorentzian full width at half maximum, in Hz.

    Raises
    ======
    ValueError
        If the input is not a finite time-domain decay or the fit does not
        return a finite positive linewidth.
    """
    if s.get_ft_prop(direct):
        raise ValueError("fit_envelope requires time-domain data")
    envelope = abs(s[direct:(0, None)]).mean_all_but([direct])
    if mult_two:
        envelope[direct, 0] *= 2
    if (
        envelope.data.size < 3
        or not np.all(np.isfinite(envelope.data))
        or envelope.data.max() <= 0
    ):
        raise ValueError("fit_envelope requires a finite nonzero decay")
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
        envelope.set_to_guess()
        orig_guess = envelope.eval()
    envelope.fit()
    new_guess = envelope.output()
    if fl:
        fl.push_marker()
        fl.next(plot_name)
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
        envelope.set_to_guess()
        points_over = (
            envelope[direct:(0, t_at_exp_end)]
            - envelope.eval()[direct:(0, t_at_exp_end)]
        )
        points_over[lambda x: x < 0] = 0
        points_over.run(lambda x: np.sqrt(abs(x) ** 2)).mean()
        amount_over[j] = points_over.data.item()
    # }}}
    lamb = "$\\lambda_L$"
    env_expansion = psp.nddata(amount_over / amount_over.max(), [-1], [lamb])
    env_expansion.setaxis(lamb, lw_range).set_units(lamb, "Hz")
    env_expansion.name("norm of points\noutside envelope")
    if fl and show_expanding_envelope:
        fl.next("expanding envelope", figsize=r_[1, 0.3] * full_width)
        fl.plot(env_expansion)
    norm_max = env_expansion.max().item().real
    norm_min = env_expansion.min().item().real
    opt_lambda = env_expansion.invinterp(
        lamb, norm_min * (1 - threshold) + norm_max * threshold, kind="linear"
    )
    logging.debug(
        "opt_lambda",
        opt_lambda,
        "at",
        norm_min * (1 - threshold) + norm_max * threshold,
        "out of",
        lw_range,
    )
    if fl and show_expanding_envelope:
        fl.plot(opt_lambda, "o")
    new_guess.update(lambda_L=opt_lambda.getaxis(lamb).item().real)
    envelope.set_guess(new_guess)
    envelope.set_to_guess()
    if fl:
        fl.next(plot_name)
    env_out = envelope.output()
    if fl:
        fl.plot(
            envelope.eval() / env_out["A"], lw=1.1, label=r"optimal envelope"
        )
    # lsq
    if fl:
        fl.plot(
            L2G(env_out["lambda_L"], criterion="energy")(s.fromaxis(direct))
        )
        fl.pop_marker()
    homogeneous_linewidth = env_out["lambda_L"]
    if (
        not np.isscalar(homogeneous_linewidth)
        or not np.isfinite(homogeneous_linewidth)
        or homogeneous_linewidth <= 0
    ):
        raise ValueError(
            "fit_envelope did not return a finite positive linewidth"
        )
    return float(homogeneous_linewidth)


def L2G(
    lambda_L,
    criterion="energy",
):
    assert np.isscalar(lambda_L)
    if criterion == "energy":
        # equal energy:
        return lambda t2: np.exp(
            -0.5 * pi**3 * lambda_L**2 * t2**2 + pi * lambda_L * abs(t2)
        )
    elif criterion == "width":
        # equal linewidth
        return lambda t2: np.exp(
            pi * lambda_L * (-pi * lambda_L * t2**2 / 4 / np.log(2) + abs(t2))
        )
