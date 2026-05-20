#!/usr/bin/env python3
"""
fit ESR to Lorentzian basis
===========================

λ_L is FWHM, not HWHM.
For Lorentzian absorption:
    L(B; B_c, λ_L) has HWHM = λ_L/2

Fourier rule:
    ℱ{L}(u) = exp(−π λ_L |u|) exp(−i 2π u B_c)

Derivative rule:
    ℱ{∂L/∂B} = i 2π u ℱ{L}

Therefore:
    A_u(u, B_c, λ_L)
      = i 2π u exp(−π λ_L |u|) exp(−i 2π u B_c)

This uses pyspecdata broadcasting over:
    u        conjugate to B₀
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
from numpy import (
    r_,
    pi,
    exp,
    logspace,
    sqrt,
    log10,
    array,
    cumsum,
    searchsorted,
    unique,
    where,
)
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
    u = d.C.ift(Bname).fromaxis(Bname)
    A = (
        -1j
        * 2
        * pi
        * u
        * exp(-pi * lambda_L * abs(u))
        * exp(1j * 2 * pi * u * center)
    )
    A[Bname,0] *= 0.5 # u=0 Heaviside issue
    A.ft(Bname)
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
    H = concat(
        [
            nddata(
                hermval(x_scaled, [0] * order + [1])
                / sqrt(2**order * factorial(order)),
                Bname,
            )
            .setaxis(Bname, x)
            .set_units(Bname, d.get_units(Bname))
            for order in range(n_hermite)
        ],
        "basis",
    )
    # Hermites use the standard H_n/sqrt(2^n n!) relative normalization,
    # then get ~5x the Lorentzian L2 norm so the L1 path prefers smooth
    # baseline terms over very wide Lorentzians.
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

    # The imaginary component is transform-roundoff; the solver boundary below
    # uses the real part.  Keep A as nddata until that boundary.

    # Collapse the physical coefficient grid only after constructing the basis.
    # The coefficient vector c still indexes individual (center, λ_L) components.
    A.smoosh(["center", "lambda_L"], "basis")
    n_lorentzian_basis = A.shape["basis"]
    lorentzian_lambda_L = A.getaxis("basis")["lambda_L"].copy()
    A.setaxis("basis", r_[0:n_lorentzian_basis]).set_units("basis", None)
    A = concat(
        [
            A,
            build_hermite_baseline_basis(
                d,
                Bname,
                baseline_norm_ratio
                * sqrt((abs(A) ** 2).sum(Bname)).data.max(),
            ).setaxis(
                "basis",
                r_[n_lorentzian_basis : n_lorentzian_basis + n_hermite],
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
    coef_lorentzian = coef_show.C["basis", 0:n_lorentzian_basis].data
    lambda_levels = unique(lorentzian_lambda_L)
    lambda_mass = array(
        [
            coef_lorentzian[lorentzian_lambda_L == this_lambda].sum()
            for this_lambda in lambda_levels
        ]
    )
    if lambda_mass.sum() > 0:
        coef_integral = cumsum(lambda_mass)
        cut1 = searchsorted(coef_integral, coef_integral[-1] / 3) + 1
        cut2 = searchsorted(coef_integral, 2 * coef_integral[-1] / 3) + 1
    else:
        cut1, cut2 = len(lambda_levels) // 3, 2 * len(lambda_levels) // 3
    cut1 = min(max(cut1, 1), len(lambda_levels) - 2)
    cut2 = min(max(cut2, cut1 + 1), len(lambda_levels) - 1)
    lambda_buckets = [
        ("narrow", lambda_levels[:cut1]),
        ("medium", lambda_levels[cut1:cut2]),
        ("broad", lambda_levels[cut2:]),
    ]
    bucket_parts = []
    for bucket_name, bucket_lambdas in lambda_buckets:
        bucket_idx = where(
            (lorentzian_lambda_L >= bucket_lambdas[0])
            & (lorentzian_lambda_L <= bucket_lambdas[-1])
        )[0]
        bucket_coef = coef_show.C["basis", 0:n_lorentzian_basis]
        bucket_coef.data[:] = 0
        bucket_coef.data[bucket_idx] = coef_lorentzian[bucket_idx]
        bucket_parts.append(
            (
                bucket_name,
                bucket_lambdas[0],
                bucket_lambdas[-1],
                A.C["basis", 0:n_lorentzian_basis].along("basis") @ bucket_coef,
            )
        )
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

    fl.next("fit at end of path")
    plot(d, label="data", alpha=0.7)
    plot(fit_show, label="full reconstruction", alpha=0.7)
    plot(d - fit_show, label="full residual", alpha=0.7)
    plot(baseline, label="baseline", alpha=0.7)
    plot(baseline_subtracted, label="baseline subtracted", alpha=0.7)
    for bucket_name, lambda_min, lambda_max, bucket_part in bucket_parts:
        plot(
            bucket_part,
            label=(
                f"{bucket_name} λ_L "
                f"{1e3 * lambda_min:.3g}-{1e3 * lambda_max:.3g} mT"
            ),
            alpha=0.2,
        )
    legend()
# }}}
