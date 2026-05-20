#!/usr/bin/env python3
"""
fit ESR to Lorentzian basis
===========================

λ_L is FWHM, not HWHM.
For Lorentzian absorption:
    L(B; B_c, λ_L) has HWHM = λ_L/2

Field-domain derivative kernel:
    A(B; B_c, λ_L)
      = real(-i 2π (1 + i (B − B_c)/λ_L)^-2)

This uses pyspecdata broadcasting over:
    B₀       field axis
    center   Lorentzian center field B_c
    λ_L      Lorentzian FWHM

The SVD compression works like this:
====================================

Original constrained problem:
    minimize ½ ‖A c − y‖₂²
    subject to c ≥ 0 and 1ᵀc ≤ τ

SVD:
    A ≈ Uᵣ Σᵣ Vᵣᵀ

Define:
    Ã = Σᵣ Vᵣᵀ
    ỹ = Uᵣᵀ y

Compressed problem:
    minimize ½ ‖Ã c − ỹ‖₂²
    subject to c ≥ 0

Positive LARS computes the LASSO path:
    minimize ½ ‖Ã c − ỹ‖₂² + α 1ᵀc
    with c ≥ 0

Reading off 1ᵀc along that path gives the desired residual-vs-L1 curve.

"""
# vim: set foldmethod=marker :

from pyspecdata import *
from numpy import r_, pi, logspace, sqrt, log10, real
from numpy.polynomial.hermite import hermval
from matplotlib.pyplot import title, xlabel, ylabel, legend
from sklearn.linear_model import lars_path
from math import factorial
import os

# {{{ changeable parameters
Bname = "$B_0$"
esr_file = "15N_S175R1a_pR_DHPC_today_200304.DSC"
# Use a tiny basis first so that we can inspect the functions and debug fast.
preview_n_center = 7
preview_n_lambda_L = 4
# The dense basis can be made larger again once we know everything is correct.
fit_n_center = 80
fit_n_lambda_L = 8
n_hermite = 4
baseline_norm_ratio = 1
# }}}


def build_lorentzian_basis(
    d,
    Bname,
    n_center,
    n_lambda_L,
    center_limits=None,
):
    x = d.getaxis(Bname)
    # go for 5x the pixel size (decayed to 0 at end), b/c otherwise, we get weird discretization issues
    lambda_L_limits = ((x[1]-x[0])*5, (x[-1]-x[0])/2)
    if center_limits is None:
        center_limits = (x[0], x[-1])
    lambda_L = nddata(
        logspace(log10(lambda_L_limits[0]), log10(lambda_L_limits[1]), n_lambda_L),
        "lambda_L",
    ).set_units("lambda_L", "T")
    center = nddata(r_[0:1 : n_center * 1j], "center")
    # Assume resonances are captured inside the acquired window: for each
    # linewidth, only generate centers at least λ_L/3 from either edge.
    # Broader off-window structure belongs in the Hermite baseline.
    center = (
        center_limits[0]
        + lambda_L / 3
        + center * (center_limits[1] - center_limits[0] - 2 * lambda_L / 3)
    )
    B = d.fromaxis(Bname) - center
    A = real(-1j * 2 * pi * (1 + 1j * B / lambda_L) ** -2)
    A.setaxis(Bname, x)
    A.set_units(Bname, d.get_units(Bname))
    # Normalize each column:
    #     ‖A_i‖₂ = 1
    #
    # Otherwise the L1 constraint would prefer some linewidths simply because
    # the column norm changes with λ_L.
    A /= sqrt((abs(A) ** 2).sum(Bname))
    return A


def build_hermite_baseline_basis(d, Bname, target_norm):
    x = d.getaxis(Bname)
    x_scaled = 2 * (x - x.mean()) / (x[-1] - x[0])
    hermites = [
        nddata(
            hermval(x_scaled, [0] * order + [1])
            / sqrt(2**order * factorial(order)),
            Bname,
        )
        .setaxis(Bname, x)
        .set_units(Bname, d.get_units(Bname))
        for order in range(n_hermite)
    ]
    H = concat(hermites + [-j for j in hermites], "basis")
    # Hermites use the standard H_n/sqrt(2^n n!) relative normalization, then
    # get boosted so the L1 path prefers smooth baseline terms over very wide
    # Lorentzians.
    H *= target_norm / sqrt((abs(H) ** 2).sum(Bname))
    return H


init_logging(level="info")

# {{{ load a real cw ESR spectrum

d = find_file(
    esr_file, exp_type="francklab_esr/Sam"
)
d.chunk_auto("harmonic", "phase")
d = d["harmonic", 0]["phase", 0]

d[Bname] *= 1e-4
d.set_units(Bname, "T")
d.set_ft_initial(Bname, "f").set_ft_prop(Bname, "time_not_aliased")
# }}}

# {{{ plots

with figlist_var() as fl:
    # {{{ construct dense Lorentzian-derivative basis using labeled broadcasting
    preview_A = build_lorentzian_basis(
        d, Bname, preview_n_center, preview_n_lambda_L
    )
    A = build_lorentzian_basis(d, Bname, fit_n_center, fit_n_lambda_L)
    print("preview basis shape", preview_A.shape)

    # Keep A as nddata until the solver boundary below.

    # Collapse the physical coefficient grid only after constructing the basis.
    # The coefficient vector c still indexes individual (center, λ_L) components.
    A.smoosh(["center", "lambda_L"], "basis")
    n_lorentzian_basis = A.shape["basis"]
    A.setaxis("basis", r_[0:n_lorentzian_basis]).set_units("basis", None)
    H = build_hermite_baseline_basis(
        d,
        Bname,
        baseline_norm_ratio * sqrt((abs(A) ** 2).sum(Bname)).data.max(),
    )
    A = concat(
        [
            A,
            H.setaxis(
                "basis",
                r_[n_lorentzian_basis : n_lorentzian_basis + H.shape["basis"]],
            ),
        ],
        "basis",
    )
    # }}}

    fl.next("reduced basis preview")
    for center_j in range(preview_n_center):
        for lambda_j in range(preview_n_lambda_L):
            plot(
                preview_A["center", center_j]["lambda_L", lambda_j],
                alpha=0.35,
                human_units=False,
            )
    title("reduced Lorentzian-derivative basis")

    # {{{ SVD-compress residual coordinates
    U, Sigma, Vh = A.C.svd(Bname, "basis")
    U.set_units(Bname, d.get_units(Bname))
    A_tilde = Sigma * Vh
    y_tilde = U.C.reorder(["SV", Bname]).along(Bname) @ d
    print(
        "during compression, d was reduced from",
        d.shape,
        "to",
        y_tilde.shape,
    )
    # }}}

    # {{{ positive LARS solver boundary
    print("beginning LARS path")
    alphas, active, coefs = lars_path(
        A_tilde.C.reorder(["SV", "basis"]).data.real,
        y_tilde.C.reorder("SV").data.real,
        method="lasso",
        positive=True,
        max_iter=400,
        return_path=True,
    )
    print("done with LARS path")
    # Immediately wrap solver output back into labeled nddata.
    coef_path = nddata(coefs, ["basis", "alpha"])
    coef_path.setaxis("basis", A_tilde.getaxis("basis"))
    coef_path.setaxis("alpha", alphas)
    # }}}

    # {{{ evaluate path with pyspecdata algebra
    # Show the least-regularized point in the path.
    coef_show = coef_path.C["alpha", -1]
    baseline = (
        A.C["basis", n_lorentzian_basis:].along("basis")
        @ coef_show.C["basis", n_lorentzian_basis:]
    )
    baseline_subtracted = (
        A.C["basis", 0:n_lorentzian_basis].along("basis")
        @ coef_show.C["basis", 0:n_lorentzian_basis]
    )
    weighted_kernel = (A * coef_show).C.reorder(["basis", Bname])
    fit_show = baseline + baseline_subtracted
    # }}}

    fl.next("positive LARS path")
    plot(
        coef_path.C.sum("basis").data,
        sqrt(
            (
                abs((A_tilde.C.along("basis") @ coef_path) - y_tilde) ** 2
            ).sum("SV")
        ).data,
        "o-",
    )
    xlabel("positive L1 mass 1ᵀc")
    ylabel("compressed residual norm ‖Ãc − ỹ‖₂")
    title("positive Lorentzian/Hermite LASSO path")

    fl.next("weighted basis functions")
    fl.image(weighted_kernel, interpolation="auto")
    title("basis functions times fitted coefficients")

    fl.next("fit at end of path")
    plot(d, label="data", alpha=0.7)
    plot(fit_show, label="full reconstruction", alpha=0.7)
    plot(d - fit_show, label="full residual", alpha=0.7)
    plot(baseline, label="baseline", alpha=0.7)
    plot(baseline_subtracted, label="baseline subtracted", alpha=0.7)
    legend()
# }}}
