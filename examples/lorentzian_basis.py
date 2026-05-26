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
from numpy import r_, pi, logspace, sqrt, log10, real, empty, ones_like
from numpy.polynomial.hermite import hermval
from matplotlib.pyplot import title, xlabel, ylabel, legend, subplot, axvline
from sklearn.linear_model import lars_path
from math import factorial, ceil
import os
import warnings

# {{{ changeable parameters
Bname = "$B_0$"
esr_file = "15N_S175R1a_pR_DHPC_today_200304.DSC"
# Use a tiny basis first so that we can inspect the functions and debug fast.
preview_n_center = 7
preview_n_lambda_L = 8
# The dense basis can be made larger again once we know everything is correct.
fit_n_center = 200
fit_n_lambda_L = 8
n_hermite = 5 
hermite_amplitude_scale = 1
lorentzian_B_range = (0.3468, 0.356)
lambda_frac_from_edge = 0  # 0 keeps all centers inside lorentzian_B_range
# Weight the stacked outside-region rows: this makes Hermite-only baseline
# mismatch outside the active spectrum cost more than ordinary full-fit RMS.
baseline_region_rms_multiplier = 5
coef_threshold_frac = 1e-2
lars_initial_max_iter = 5000
lars_max_iter_limit = 50000
lars_residual_target_multiplier = 2
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
    active_width = center_limits[1] - center_limits[0]
    if lambda_frac_from_edge > 0:
        lambda_L_max = active_width / (2 * lambda_frac_from_edge)
    else:
        lambda_L_max = active_width
    if lambda_L_max <= lambda_L_min:
        raise ValueError(
            "No Lorentzian linewidths fit inside lorentzian_B_range with "
            f"lambda_frac_from_edge={lambda_frac_from_edge}"
        )
    lambda_L_limits = (
        lambda_L_min,
        lambda_L_max,
    )
    lambda_L_values = logspace(
        log10(lambda_L_limits[0]), log10(lambda_L_limits[1]), n_lambda_L
    )
    center_values = r_[center_limits[0] : center_limits[1] : n_center * 1j]
    center_pairs = center_values[:, None] + 0 * lambda_L_values[None, :]
    lambda_L_pairs = 0 * center_values[:, None] + lambda_L_values[None, :]
    # Use one fixed center grid for all linewidths.
    # Optionally reject centers close to the active-region edge, but
    # default to no rejection now that the constant offset is removed
    # explicitly below.
    if lambda_frac_from_edge > 0:
        valid_pairs = (
            center_pairs
            >= center_limits[0] + lambda_frac_from_edge * lambda_L_pairs
        ) & (
            center_pairs
            <= center_limits[1] - lambda_frac_from_edge * lambda_L_pairs
        )
    else:
        valid_pairs = ones_like(center_pairs, dtype=bool)
    if not valid_pairs.any():
        raise ValueError("No Lorentzian center/linewidth pairs survived filtering")
    basis_axis = empty(
        valid_pairs.sum(), dtype=[("center", float), ("lambda_L", float)]
    )
    basis_axis["center"] = center_pairs[valid_pairs]
    basis_axis["lambda_L"] = lambda_L_pairs[valid_pairs]
    center = nddata(basis_axis["center"], "basis").setaxis("basis", basis_axis)
    lambda_L = nddata(basis_axis["lambda_L"], "basis").setaxis(
        "basis", basis_axis
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


def mark_active_region(active_B_range, ax=None):
    if active_B_range is None:
        return
    for edge in active_B_range:
        if ax is None:
            axvline(edge, color="k", linestyle=(0, (1, 2)), linewidth=1)
        else:
            ax.axvline(edge, color="k", linestyle=(0, (1, 2)), linewidth=1)


def fit_lars_path(
    fl,
    label,
    A,
    d,
    Bname,
    first_baseline_basis,
    coef_threshold_frac,
    active_B_range=None,
    baseline_region_rms_multiplier=1,
):
    basis_axis = A.getaxis("basis")
    lorentzian_mask = basis_axis < first_baseline_basis
    baseline_mask = ~lorentzian_mask
    if active_B_range is not None:
        x = d.getaxis(Bname)
        outside_active = (x < active_B_range[0]) | (x > active_B_range[1])
        if outside_active.any():
            A_outside = A.C[Bname, outside_active]
            # In the duplicated outside-region rows, only the Hermite baseline is
            # allowed to explain the data.  The full, original rows still measure
            # the overall reconstruction RMS.
            A_outside["basis", lorentzian_mask] *= 0
            A_outside *= baseline_region_rms_multiplier
            d_outside = d.C[Bname, outside_active] * baseline_region_rms_multiplier
            A_fit = concat([A, A_outside], Bname)
            d_fit = concat([d, d_outside], Bname)
            print(
                f"{label}: added",
                outside_active.sum(),
                "outside-active-region baseline rows with multiplier",
                baseline_region_rms_multiplier,
            )
        else:
            A_fit = A
            d_fit = d
            print(f"{label}: no outside-active-region rows to add")
    else:
        A_fit = A
        d_fit = d
    # {{{ SVD-compress residual coordinates
    U, Sigma, Vh = A_fit.C.svd(Bname, "basis")
    U.set_units(Bname, d_fit.get_units(Bname))
    A_tilde = Sigma * Vh
    y_tilde = U.C.reorder(["SV", Bname]).along(Bname) @ d_fit
    print(
        f"{label}: during compression, d was reduced from",
        d_fit.shape,
        "to",
        y_tilde.shape,
    )
    # }}}

    # {{{ positive LARS solver boundary
    X = A_tilde.C.reorder(["SV", "basis"]).data.real
    y = y_tilde.C.reorder("SV").data.real
    x = d.getaxis(Bname)
    if active_B_range is None:
        inside_active = ones_like(x, dtype=bool)
    else:
        inside_active = (x >= active_B_range[0]) & (x <= active_B_range[1])
    outside_active = ~inside_active

    def coefficient_residual_stats(coef_values):
        coef_candidate = nddata(coef_values, "basis")
        coef_candidate.setaxis("basis", A.getaxis("basis"))
        baseline_candidate = (
            A.C["basis", baseline_mask].along("basis")
            @ coef_candidate.C["basis", baseline_mask]
        )
        baseline_subtracted_candidate = (
            A.C["basis", lorentzian_mask].along("basis")
            @ coef_candidate.C["basis", lorentzian_mask]
        )
        for j in [baseline_candidate, baseline_subtracted_candidate]:
            j.data = real(j.data)
        fit_candidate = baseline_candidate + baseline_subtracted_candidate
        active_residual = (d - fit_candidate).C[Bname, inside_active]
        active_rms = sqrt((abs(active_residual) ** 2).mean(Bname)).data.item()
        if outside_active.any():
            baseline_residual = (
                d.C[Bname, outside_active]
                - baseline_candidate.C[Bname, outside_active]
            )
            baseline_rms = sqrt(
                (abs(baseline_residual) ** 2).mean(Bname)
            ).data.item()
        else:
            baseline_rms = active_rms / lars_residual_target_multiplier
        return {
            "active_rms": active_rms,
            "baseline_rms": baseline_rms,
            "target_rms": lars_residual_target_multiplier * baseline_rms,
        }

    lars_max_iter = max(lars_initial_max_iter, 3 * X.shape[1])
    while True:
        print(f"{label}: beginning LARS path with max_iter", lars_max_iter)
        alphas, active, coefs = lars_path(
            X,
            y,
            method="lasso",
            positive=True,
            max_iter=lars_max_iter,
            alpha_min=0,
            return_path=True,
        )
        print(f"{label}: done with LARS path")
        final_stats = coefficient_residual_stats(coefs[:, -1])
        if final_stats["active_rms"] <= final_stats["target_rms"]:
            break
        if lars_max_iter >= lars_max_iter_limit:
            warnings.warn(
                f"{label}: LARS reached max_iter={lars_max_iter} but active "
                f"RMS {final_stats['active_rms']:.3g} remained above target "
                f"{final_stats['target_rms']:.3g}; using best available point"
            )
            break
        target_rms = max(final_stats["target_rms"], 1e-15)
        new_max_iter = min(
            lars_max_iter_limit,
            ceil(
                lars_max_iter
                * max(2, 1.5 * final_stats["active_rms"] / target_rms)
            ),
        )
        warnings.warn(
            f"{label}: final LARS active RMS {final_stats['active_rms']:.3g} "
            f"exceeded target {final_stats['target_rms']:.3g}; rerunning with "
            f"max_iter={new_max_iter}"
        )
        lars_max_iter = new_max_iter
    # LARS coefficients are still on A's original basis axis; only the
    # residual coordinates have been compressed.
    coef_path = nddata(coefs, ["basis", "alpha"])
    coef_path.setaxis("basis", A_tilde.getaxis("basis"))
    coef_path.setaxis("alpha", alphas)

    if final_stats["active_rms"] <= final_stats["target_rms"]:
        selected_idx = coefs.shape[1] - 1
        selected_stats = final_stats
        for j in range(coefs.shape[1] - 2, -1, -1):
            this_stats = coefficient_residual_stats(coefs[:, j])
            if this_stats["active_rms"] <= this_stats["target_rms"]:
                selected_idx = j
                selected_stats = this_stats
            else:
                break
    else:
        selected_idx = coefs.shape[1] - 1
        selected_stats = final_stats
        for j in range(coefs.shape[1] - 2, -1, -1):
            this_stats = coefficient_residual_stats(coefs[:, j])
            if this_stats["active_rms"] < selected_stats["active_rms"]:
                selected_idx = j
                selected_stats = this_stats
    print(
        f"{label}: selected LARS index",
        selected_idx,
        "alpha",
        alphas[selected_idx],
        "active RMS",
        selected_stats["active_rms"],
        "baseline RMS",
        selected_stats["baseline_rms"],
        "target RMS",
        selected_stats["target_rms"],
    )
    # }}}

    # {{{ evaluate path with pyspecdata algebra
    coef_show = coef_path.C["alpha", selected_idx]
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
    selected_l1 = coef_show.C.sum("basis").data.item()
    axvline(selected_l1, linestyle="--", color="k")
    plot(
        [selected_l1],
        [residual_path[selected_idx]],
        "kx",
        label="selected residual target",
    )
    xlabel(r"$\mathbf{1}^{\mathsf{T}}\mathbf{c}$")
    ylabel(
        r"$\left\|\widetilde{\mathbf{A}}\mathbf{c}"
        r"-\widetilde{\mathbf{y}}\right\|_2$"
    )
    title("positive Lorentzian/Hermite LASSO path")
    legend()

    fl.next(f"{label}: weighted basis functions")
    print(weighted_kernel.data.dtype)
    ax = subplot(1, 2, 1)
    fl.image(weighted_kernel, interpolation="auto", ax=ax, human_units=False)
    mark_active_region(active_B_range, ax=ax)
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
    mark_active_region(active_B_range)
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
# Remove the constant offset explicitly; otherwise the positive Lorentzian
# and baseline basis spend most of their effort reproducing DC.
dc_offset = d.C.mean(Bname)
print("subtracting constant offset", dc_offset.data.item())
d -= dc_offset
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

    n_lorentzian_basis = A.shape["basis"]
    lorentzian_basis_axis = A.getaxis("basis")
    lambda_L_max = lorentzian_basis_axis["lambda_L"].max()
    A.setaxis("basis", r_[0:n_lorentzian_basis]).set_units("basis", None)

    # {{{ build Hermite baseline basis: generate ± Hermite polynomial columns
    # Hermites use the standard H_n/sqrt(2^n n!) relative normalization.
    # Larger basis-function amplitude is cheaper in LARS because the same
    # signal needs less coefficient.  At hermite_amplitude_scale=1, scale the
    # Hermites so their largest absolute excursion is 1/lambda_L_max, making
    # them approximately competitive with the broadest allowed Lorentzians.
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
    hermite_target_amp = hermite_amplitude_scale / lambda_L_max
    hermite_scale = hermite_target_amp / abs(H).data.max()
    print("lambda_L max", lambda_L_max)
    print("Hermite target amplitude", hermite_target_amp)
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
    for basis_j in range(preview_A.shape["basis"]):
        plot(
            preview_A["basis", basis_j],
            alpha=0.35,
            human_units=False,
        )
    mark_active_region(lorentzian_B_range)
    title("reduced Lorentzian-derivative basis")

    full_fit = fit_lars_path(
        fl,
        "full basis",
        A,
        d,
        Bname,
        n_lorentzian_basis,
        coef_threshold_frac,
        lorentzian_B_range,
        baseline_region_rms_multiplier,
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
        lorentzian_B_range,
        baseline_region_rms_multiplier,
    )
# }}}
