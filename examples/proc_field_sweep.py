"""
hydroxyTEMPO ODNP field sweep: DNP spectrum (epsilon = 1 - E)
==============================================================

``py proc_field_sweep.py ON_FILE OFF_FILE EXP_TYPE [ON_NODE OFF_NODE]``

Loads the MW-on and MW-off field-sweep files/nodes with
``find_file(..., lookup=lookup_table)`` (dispatching on the saved
``field_sweep_v5`` postproc, so each node returns conjugated,
``tau``-shifted, FT'd along ``t2`` and along the phase cycle, with
``coherence_pathway`` / ``acq_params`` / ``field_readback_G`` as
properties), phases each node, integrates the on-resonance NMR peak for
every field, and plots ``epsilon = 1 - E`` (``E = S_on / S_off``) with a
three-line pseudo-Voigt fit via ``lmfitdata``, alongside the analytic
derivative of the fit.

Phasing (before any integration or plotting)
---------------------------------------------
Per repo convention (see ``examples/Hermitian_Phasing.py``), each node gets:

Opt. **Equal-energy apodization**
1. **zeroth-order phasing** (``zeroth_order_ph``) -- a single global phase
   for the whole node, found from the moment of inertia of the
   coherence-selected signal in the complex plane.
2. **Hermitian (first-order) phasing** (``hermitian_function_test``) -- a
   single global timing correction that centers the echo, found from the
   nu_offset-averaged, coherence-selected FID/echo.
3. **Correlation alignment**

Sign convention
----------------
``E = S_on / S_off`` is negative under the DNP lines (1H Overhauser
enhancement inverts the water-proton signal).  We plot the more
conventional ``epsilon = 1 - E`` instead, which is positive under the DNP
lines and zero off-resonance.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp
from matplotlib.gridspec import GridSpec
from pyspecProcScripts import (
    L2G,
    correl_align,
    hermitian_function_test,
    select_pathway,
    zeroth_order_ph,
)

# {{{ changeable parameters

on_file, off_file, exp_type = (
    "260720_hydroxytempo_field_sweep.h5",
    "260720_hydroxytempo_field_sweep.h5",
    "b27/field_dependent",
)
on_node, off_node = "field_sweep_1", "field_sweep_2"  # MW on / MW off
signal_pathway = {"ph1": 1}
peak_lower_thresh = 0.1
apodization = True  # Apply to both or none. Recommended.
apod_width_Hz = 1e3  # After the optimal value increasing shouldn't change
# the enhancement. Low values give weird results.
# }}}}


def phase_node_and_align(filename, nodename, apod=False, fl=None):
    """Load one node and apply zeroth-order phasing and Hermitian timing
    correction.  Optionally apply equal-energy apodization and correlation
    alignment.  Returns the phased nddata, acq_params, the center field, and
    the zeroth-order phase actually used (so the caller can share it with the
    other node)."""
    s = psd.find_file(
        filename, exp_type=exp_type, expno=nodename, lookup=prscr.lookup_table
    )
    acq_params = s.get_prop("acq_params")
    field_G = s.get_prop("field_readback_G")["NUMPY_DATA"]
    center_field_G = (
        acq_params["carrierFreq_MHz"] / acq_params["gamma_eff_MHz_G"]
    )
    s["indirect"] = (
        (field_G - center_field_G)
        * acq_params["uw_dip_center_GHz"]
        * 1e3
        / center_field_G
    )
    s.rename("indirect", "nu_offset").set_units("nu_offset", "MHz")
    s.sort("nu_offset")
    s.set_prop("coherence_pathway", signal_pathway)

    if apod:
        # {{{ Equal-energy apodization
        fig = fl.next(f"{nodename} apodized signal")
        fl.title = f"{nodename} apodized signal"
        gs = GridSpec(1, 3, figure=fig, left=0.05, right=0.95)
        psd.DCCT(
            s,
            fig,
            title="Raw Data",
            bbox=gs[0, 0],
        )
        s.ift("t2")
        s *= L2G(apod_width_Hz, criterion="energy")(s.fromaxis("t2"))
        s.ft("t2")
        psd.DCCT(
            s,
            fig,
            title="Equal Energy Apodization",
            bbox=gs[0, 1],
        )
    signal = select_pathway(s, signal_pathway)
    nu_axis = signal.getaxis("nu_offset")
    frq_center, frq_half = prscr.find_peakrange(
        signal["nu_offset" : nu_axis[abs(nu_axis).argmin()]],
        peak_lower_thresh=peak_lower_thresh,
    )
    frq_half = abs(frq_half)

    s.ift("t2")
    zeroth_phase = zeroth_order_ph(select_pathway(s, signal_pathway))
    s /= zeroth_phase
    s["t2"] -= s.getaxis("t2")[0]  # needed for hermitian_function_test
    best_shift = hermitian_function_test(
        select_pathway(s.C.mean("nu_offset"), signal_pathway)
    )
    s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
    s.ft("t2")

    def frq_mask(this_s):
        repeat_dim = (
            "repeats" if "repeats" in this_s.dimlabels else "nu_offset"
        )
        nu_center = (
            select_pathway(this_s, signal_pathway)
            .mean(repeat_dim)
            .argmax("t2")
        )
        return this_s * np.exp(
            -((this_s.fromaxis("t2") - nu_center) ** 2) / (4 * frq_half**2)
        )

    def coherence_unmask_fn(coh_array):
        thisslice = coh_array
        for j, (k, v) in enumerate(signal_pathway.items()):
            if j == len(signal_pathway) - 1:
                thisslice[k, v] = 1
            else:
                thisslice = thisslice[k, v]
        return coh_array

    repeat_sign = select_pathway(s, signal_pathway).C.real.sum("t2")
    repeat_sign = repeat_sign.run(np.sign)
    opt_shift = correl_align(
        s * repeat_sign,
        frq_mask_fn=frq_mask,
        coherence_unmask_fn=coherence_unmask_fn,
        repeat_dims="nu_offset",
        max_shift=frq_half,
    )
    s.ift("t2").ift(list(signal_pathway))
    s *= np.exp(-1j * 2 * np.pi * opt_shift * s.fromaxis("t2"))
    s.ft("t2").ft(list(signal_pathway.keys()))

    if apod:
        psd.DCCT(
            s,
            fig,
            title="After apodization\nand alignment",
            bbox=gs[0, 2],
        )

    return s, acq_params, center_field_G, frq_center, frq_half


with psd.figlist_var() as fl:
    fl.basename = on_file

    (
        off,
        acq_params,
        center_field_G,
        off_frq_center,
        off_frq_half,
    ) = phase_node_and_align(off_file, off_node, apod=apodization, fl=fl)
    on, _, _, on_frq_center, on_frq_half = phase_node_and_align(
        on_file, on_node, apod=apodization, fl=fl
    )

    on_band = select_pathway(on, signal_pathway)[
        "t2" : (on_frq_center - on_frq_half, on_frq_center + on_frq_half)
    ].integrate("t2")
    off_band = select_pathway(off, signal_pathway)[
        "t2" : (off_frq_center - off_frq_half, off_frq_center + off_frq_half)
    ].integrate("t2")

    # {{{ epsilon = 1 - E, E = S_on / S_off, off node interpolated onto the
    #     on-node field axis.
    nu_on = np.asarray(on_band["nu_offset"], dtype=float)
    nu_off = np.asarray(off_band["nu_offset"], dtype=float)
    off_order = np.argsort(nu_off)
    off_on_axis = np.interp(
        nu_on, nu_off[off_order], off_band.data.real[off_order]
    )
    thermal_level = np.median(off_on_axis)
    # DNP data is negative wrt thermals
    if np.median(on_band.data.real) * thermal_level > 0:
        on_band *= -1
    epsilon = on_band.C
    epsilon.data = 1 - on_band.data.real / thermal_level
    epsilon.set_error(None)  # None for now. In the future we can add errors
    # from the goodness of the fit.

    epsilon.name("epsilon")
    # }}}

    # {{{ three-line pseudo-Voigt fit (shared center + 14N hyperfine spacing)
    nu_offset, E_0, sigma, lambda_L, eta, nu_c, A_hf, A_0, A_1, A_2 = (
        sp.symbols(
            "nu_offset E_0 sigma lambda_L eta nu_c A_hf A_0 A_1 A_2",
            real=True,
        )
    )
    fitdata = psd.lmfitdata(epsilon)
    fitdata.functional_form = E_0 + sum(
        A
        * (
            eta * lambda_L**2 / ((nu_offset - c) ** 2 + lambda_L**2)
            + (1 - eta) * sp.exp(-((nu_offset - c) ** 2) / (2 * sigma**2))
        )
        for A, c in (
            (A_0, nu_c - A_hf),
            (A_1, nu_c),
            (A_2, nu_c + A_hf),
        )
    )
    depth = float(abs(epsilon.data).max())
    sweep = float(np.ptp(nu_on))
    fitdata.set_guess(
        E_0=dict(value=0.0, min=-depth, max=depth),
        nu_c=dict(value=0.0, min=-10.0, max=10.0),
        A_hf=dict(value=45.0, min=30.0, max=60.0),  # 14N hyperfine spacing
        A_0=dict(value=depth, min=-20 * depth, max=20 * depth),
        A_1=dict(value=depth, min=-20 * depth, max=20 * depth),
        A_2=dict(value=depth, min=-20 * depth, max=20 * depth),
        sigma=dict(value=sweep / 12, min=0.5, max=sweep / 2),
        lambda_L=dict(value=sweep / 12, min=0.5, max=sweep / 2),
        eta=dict(value=0.5, min=0.01, max=0.99),
    )
    # stage 1: line positions frozen; stage 2: released from stage-1 result
    for guess_name, guess_par in fitdata.guess_parameters.items():
        guess_par.vary = guess_name not in ("nu_c", "A_hf")
    fitdata.fit(use_jacobian=False)
    fitdata.guess_parameters = fitdata.fit_parameters
    for guess_par in fitdata.guess_parameters.values():
        guess_par.vary = True
    fitdata.fit(use_jacobian=False)
    # }}}

    # {{{ analytic derivative of the fitted curve -- comparable to a
    #     conventional (derivative-mode) EPR spectrum: Eq. S7-S9 in the
    #     supplement write the DNP transition rates in terms of the EPR
    #     *absorption* lineshape phi(nu), so our fitted epsilon curve is
    #     itself absorption-shaped, and its derivative is what a field-
    #     modulated CW-EPR spectrometer would display-like.  We differentiate
    #     the noise-free fitted expression (substituting fitted parameter
    #     values via their sympy symbols, not the string names output()
    #     returns) rather than the raw data, since numerically
    #     differentiating scattered points would just amplify noise.
    symbol_table = {
        "E_0": E_0,
        "sigma": sigma,
        "lambda_L": lambda_L,
        "eta": eta,
        "nu_c": nu_c,
        "A_hf": A_hf,
        "A_0": A_0,
        "A_1": A_1,
        "A_2": A_2,
    }
    fitted_expr = fitdata.functional_form.subs(
        {symbol_table[k]: v for k, v in fitdata.output().items()}
    )
    derivative_func = sp.lambdify(
        nu_offset, sp.diff(fitted_expr, nu_offset), "numpy"
    )
    nu_fine = np.linspace(nu_on.min(), nu_on.max(), 500)
    derivative = psd.nddata(derivative_func(nu_fine), ["nu_offset"])
    derivative.setaxis("nu_offset", nu_fine).set_units("nu_offset", "MHz")
    derivative.name("d(epsilon)/d(nu_offset)")
    # }}}

    fig, (ax_eps, ax_deriv) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))
    fl.next("hydroxyTEMPO ODNP field sweep", fig=fig)
    psd.plot(
        epsilon,
        "o",
        ax=ax_eps,
        color="k",
        markersize=3,
        label="DNP curve",
    )
    psd.plot(
        fitdata.eval(500),
        ax=ax_eps,
        color="r",
        alpha=0.8,
        label="3-line pseudo-Voigt",
    )
    nu_c_fit, A_hf_fit = fitdata.output("nu_c"), fitdata.output("A_hf")
    line_centers = (nu_c_fit - A_hf_fit, nu_c_fit, nu_c_fit + A_hf_fit)
    epr_shift_per_G = acq_params["uw_dip_center_GHz"] * 1e3 / center_field_G
    middle_peak_nmr_MHz = (
        center_field_G + nu_c_fit / epr_shift_per_G
    ) * acq_params["gamma_eff_MHz_G"]
    for this_center in line_centers:
        ax_eps.axvline(this_center, ls=":", color="grey", alpha=0.5)
    ax_eps.axhline(0, lw=0.6, color="k")
    ax_eps.set_ylabel(r"$\epsilon  = 1 - E$")
    ax_eps.legend()
    ax_eps.grid(alpha=0.3)

    psd.plot(
        derivative,
        ax=ax_deriv,
        color="b",
        label=r"$d\epsilon/d\nu$",
    )
    for this_center in line_centers:
        ax_deriv.axvline(this_center, ls=":", color="grey", alpha=0.5)
    ax_deriv.axhline(0, lw=0.6, color="k")
    ax_deriv.set_xlabel(r"EPR offset $\nu-\nu_\mathrm{center}$ / MHz")
    ax_deriv.set_ylabel("First Derivative Spectrum")
    ax_deriv.legend()
    ax_deriv.grid(alpha=0.3)
    fig.tight_layout()

    print(f"center field: {center_field_G:#0.8g} G")
    print(f"peak epsilon: {epsilon.data.max():#0.4g}")
    print(
        f"14N hyperfine spacing: {A_hf_fit:#0.4g} MHz;  centers (MHz): "
        f"{nu_c_fit - A_hf_fit:#0.4g}, {nu_c_fit:#0.4g}, "
        f"{nu_c_fit + A_hf_fit:#0.4g}"
    )
    print(
        "NMR/EPR ratio: "
        f"{middle_peak_nmr_MHz / acq_params['uw_dip_center_GHz']:#0.8g}"
        " MHz/GHz "
    )
