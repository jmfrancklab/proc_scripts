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
from matplotlib.pyplot import title, xlabel, ylabel, legend, subplot, axvline
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
n_hermite = 2
baseline_cost_multiplier = 10
lorentzian_B_range = (0.344, 0.358)
lambda_frac_from_edge = 5  # prevent lopsided contributions
coef_threshold_frac = 1e-2
# }}}


def build_lorentzian_basis(
    d,
    Bname,
    n_center,
    n_lambda_L,
    center_limits=None,
):
    x = d.getaxis(Bname)
    if center_limits is None:
        center_limits = (x[0], x[-1])
    # go for 5x the pixel size (decayed to 0 at end), b/c otherwise, we get weird discretization issues
    lambda_L_min = (x[1] - x[0]) * 5
    lambda_L_max = (center_limits[1] - center_limits[0]) / (
        2 * lambda_frac_from_edge
    )
    if lambda_L_max <= lambda_L_min:
        raise ValueError(
            "No Lorentzian linewidths fit inside lorentzian_B_range with "
            f"lambda_frac_from_edge={lambda_frac_from_edge}"
        )
    lambda_L_limits = (
        lambda_L_min,
        lambda_L_max,
    )
    lambda_L = nddata(
        logspace(
            log10(lambda_L_limits[0]), log10(lambda_L_limits[1]), n_lambda_L
        ),
        "lambda_L",
    ).set_units("lambda_L", "T")
    center = nddata(r_[0 : 1 : n_center * 1j], "center")
    # Assume resonances are captured inside the acquired window: for each
    # linewidth, only generate centers at least λ_L*lambda_frac_from_edge from either edge.
    # Broader off-window structure belongs in the Hermite baseline.
    center = (
        center_limits[0]
        + lambda_L * lambda_frac_from_edge
        + center
        * (
            center_limits[1]
            - center_limits[0]
            - 2 * lambda_L * lambda_frac_from_edge
        )
    )
    # The raw derivative shape has linewidth-independent peak amplitude.
    # Physically, equal ESR spin count would give peak-to-peak amplitude
    # proportional to 1/lambda_L**2.  Multiplying that physical scale by
    # lambda_L intentionally favors one broad Lorentzian over a collection of
    # narrow Lorentzians, so the implemented compromise divides by lambda_L.
    A = (
        real(-1j * (1 + 1j * (d.fromaxis(Bname) - center) / lambda_L) ** -2)
        / lambda_L
    )
    A.setaxis(Bname, x)
    A.set_units(Bname, d.get_units(Bname))
    return A


def fit_lars_path(
    fl,
    label,
    A,
    d,
    Bname,
    first_baseline_basis,
    coef_threshold_frac,
):
    # {{{ SVD-compress residual coordinates
    U, Sigma, Vh = A.C.svd(Bname, "basis")
    U.set_units(Bname, d.get_units(Bname))
    A_tilde = Sigma * Vh
    y_tilde = U.C.reorder(["SV", Bname]).along(Bname) @ d
    print(
        f"{label}: during compression, d was reduced from",
        d.shape,
        "to",
        y_tilde.shape,
    )
    # }}}

    # {{{ positive LARS solver boundary
    print(f"{label}: beginning LARS path")
    alphas, active, coefs = lars_path(
        A_tilde.C.reorder(["SV", "basis"]).data.real,
        y_tilde.C.reorder("SV").data.real,
        method="lasso",
        positive=True,
        max_iter=400,
        return_path=True,
    )
    print(f"{label}: done with LARS path")
    # LARS coefficients are still on A's original basis axis; only the
    # residual coordinates have been compressed.
    coef_path = nddata(coefs, ["basis", "alpha"])
    coef_path.setaxis("basis", A_tilde.getaxis("basis"))
    coef_path.setaxis("alpha", alphas)
    # }}}

    # {{{ evaluate path with pyspecdata algebra
    coef_show = coef_path.C["alpha", -1]
    basis_axis = A.getaxis("basis")
    lorentzian_mask = basis_axis < first_baseline_basis
    baseline_mask = ~lorentzian_mask
    baseline = (
        A.C["basis", baseline_mask].along("basis")
        @ coef_show.C["basis", baseline_mask]
    )
    baseline_subtracted = (
        A.C["basis", lorentzian_mask].along("basis")
        @ coef_show.C["basis", lorentzian_mask]
    )
    weighted_kernel = (A * coef_show).C.reorder(["basis", Bname])
    # SVD/dot bookkeeping can leave complex dtype; the original field-domain
    # kernel and ESR data are real.
    for j in [baseline, baseline_subtracted, weighted_kernel]:
        j.data = real(j.data)
    fit_show = baseline + baseline_subtracted
    coef_amplitudes = abs(coef_show).data
    coef_cutoff = coef_threshold_frac * coef_amplitudes[lorentzian_mask].max()
    l1_path = coef_path.C.sum("basis").data
    residual_path = sqrt(
        (abs((A_tilde.C.along("basis") @ coef_path) - y_tilde) ** 2).sum("SV")
    ).data
    # }}}

    fl.next(f"{label}: positive LARS path")
    plot(l1_path, residual_path, "o-")
    axvline(coef_show.C.sum("basis").data.item(), linestyle="--", color="k")
    xlabel(r"$\mathbf{1}^{\mathsf{T}}\mathbf{c}$")
    ylabel(
        r"$\left\|\widetilde{\mathbf{A}}\mathbf{c}"
        r"-\widetilde{\mathbf{y}}\right\|_2$"
    )
    title("positive Lorentzian/Hermite LASSO path")

    fl.next(f"{label}: weighted basis functions")
    print(weighted_kernel.data.dtype)
    ax = subplot(1, 2, 1)
    fl.image(weighted_kernel, interpolation="auto", ax=ax)
    ax.set_title("basis functions times fitted coefficients")
    ax = subplot(1, 2, 2)
    ax.semilogy(basis_axis, coef_amplitudes, ".")
    ax.axhline(coef_cutoff, linestyle="--", color="k")
    ax.set_xlabel("original basis index")
    ax.set_ylabel("coefficient amplitude")
    ax.set_title("fit coefficient amplitudes")

    fl.next(f"{label}: fit at end of path")
    plot(d, label="data", alpha=0.7)
    plot(fit_show, label="full reconstruction", alpha=0.7)
    plot(d - fit_show, label="full residual", alpha=0.7)
    plot(baseline, label="baseline", alpha=0.7)
    plot(baseline_subtracted, label="baseline subtracted", alpha=0.7)
    legend()
    return {
        "coef_show": coef_show,
        "coef_cutoff": coef_cutoff,
        "lorentzian_mask": lorentzian_mask,
    }


init_logging(level="info")

# {{{ load a real cw ESR spectrum

d = find_file(esr_file, exp_type="francklab_esr/Sam")
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
        d, Bname, preview_n_center, preview_n_lambda_L, lorentzian_B_range
    )
    A = build_lorentzian_basis(
        d, Bname, fit_n_center, fit_n_lambda_L, lorentzian_B_range
    )
    print("note that this is real, as it should be", A.data.dtype)
    print("preview basis shape", preview_A.shape)

    # Keep A as nddata until the solver boundary below.

    # Collapse the physical coefficient grid only after constructing the basis.
    # The coefficient vector c still indexes individual (center, λ_L) components.
    A.smoosh(["center", "lambda_L"], "basis")
    n_lorentzian_basis = A.shape["basis"]
    lorentzian_basis_axis = A.getaxis("basis")
    broadest_lorentzian_amp = abs(
        A[
            "basis",
            lorentzian_basis_axis["lambda_L"]
            == lorentzian_basis_axis["lambda_L"].max(),
        ]
    ).data.max()
    A.setaxis("basis", r_[0:n_lorentzian_basis]).set_units("basis", None)

    # {{{ build Hermite baseline basis: generate ± Hermite polynomial columns
    # Hermites use the standard H_n/sqrt(2^n n!) relative normalization.
    # They do not have a spin-count meaning, so price them against the broadest
    # smooth Lorentzian they might replace: scale their largest absolute
    # excursion to broadest_lorentzian_amp / baseline_cost_multiplier.  This
    # uses max amplitude rather than L2 norm because the Lorentzian expense was
    # chosen from physical peak-to-peak scaling, not equal-energy atoms.
    # Both signs are included because the solver coefficients are constrained
    # positive.
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
    hermite_scale = broadest_lorentzian_amp / (
        baseline_cost_multiplier * abs(H).data.max()
    )
    print("broadest Lorentzian amplitude", broadest_lorentzian_amp)
    print("Hermite scale factor", hermite_scale)
    H *= hermite_scale
    # }}}

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

    full_fit = fit_lars_path(
        fl,
        "full basis",
        A,
        d,
        Bname,
        n_lorentzian_basis,
        coef_threshold_frac,
    )
    basis_axis = A.getaxis("basis")
    coef_amplitudes = abs(full_fit["coef_show"]).data
    reduced_basis_mask = basis_axis >= n_lorentzian_basis
    reduced_basis_mask |= (basis_axis < n_lorentzian_basis) & (
        coef_amplitudes >= full_fit["coef_cutoff"]
    )
    print(
        "reduced basis keeps",
        reduced_basis_mask[:n_lorentzian_basis].sum(),
        "of",
        n_lorentzian_basis,
        "Lorentzian columns",
    )
    fit_lars_path(
        fl,
        "thresholded Lorentzian basis",
        A["basis", reduced_basis_mask],
        d,
        Bname,
        n_lorentzian_basis,
        coef_threshold_frac,
    )
# }}}
