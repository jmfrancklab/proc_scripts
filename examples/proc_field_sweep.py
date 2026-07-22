"""
hydroxyTEMPO ODNP field sweep: signed DNP spectrum
==================================================

``py proc_field_sweep_dnp_curve.py FILENAME EXP_TYPE``

Loads the two field-sweep nodes with ``find_file(..., lookup=lookup_table)``
(dispatching on the saved ``field_sweep_v5`` postproc, so each node returns
conjugated, ``tau``-shifted, FT'd along ``t2`` and along the phase cycle, with
``coherence_pathway`` / ``acq_params`` / ``field_readback_G`` as properties),
integrates the on-resonance NMR peak for every field, and plots the **signed**
DNP enhancement spectrum with a three-line pseudo-Voigt fit via ``lmfitdata``.

Signed processing
-----------------
``determine_sign`` assigns a per-point +-1 that scatters on the low-SNR flanks
of this long sweep, where the receiver phase drifts field-to-field.  The drift
is smooth, though, and the MW-off (thermal) node carries that same instrumental
phase with no DNP sign flips, so we build a magnitude-weighted, field-smoothed
phasor from the thermal node and rotate the MW-on band integrals into that
frame.  The real part is the signed enhancement -- negative under the DNP
lines, because 1H Overhauser enhancement inverts the water-proton signal --
with no per-point sign scatter.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp
from pyspecProcScripts import select_pathway

# {{{ changeable parameters
if len(sys.argv) == 3:
    thisfile, exp_type = sys.argv[1], sys.argv[2]
else:
    thisfile, exp_type = (
        "260720_hydroxytempo_field_sweep.h5",
        "b27/field_dependent",
    )
on_node, off_node = "field_sweep_1", "field_sweep_2"  # MW on / MW off
signal_band_Hz = 3e3  # +-band around the on-resonance NMR peak
phase_smooth_MHz = 10.0  # thermal-phase smoothing kernel width
# }}}


def band_integral(nodename):
    """Load one node, put the EPR offset (MHz) on the indirect axis, select the
    coherence pathway, and return the complex on-resonance band integral
    ``I(nu_offset)`` together with ``acq_params`` and the center field."""
    s = psd.find_file(
        thisfile, exp_type=exp_type, expno=nodename, lookup=prscr.lookup_table
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
    return (
        select_pathway(s, s.get_prop("coherence_pathway"))[
            "t2" : (-signal_band_Hz, signal_band_Hz)
        ].integrate("t2"),
        acq_params,
        center_field_G,
    )


on, acq_params, center_field_G = band_integral(on_node)
off, _, _ = band_integral(off_node)

# {{{ thermal-frame phasing -> signed enhancement E = S_on / S_off
nu = np.asarray(on["nu_offset"], dtype=float)
nu_off = np.asarray(off["nu_offset"], dtype=float)
thermal_phase = np.array(
    [
        np.angle(
            np.sum(
                np.exp(-0.5 * ((nu_off - this_nu) / phase_smooth_MHz) ** 2)
                * off.data
            )
        )
        for this_nu in nu
    ]
)
enhancement = on.C
enhancement.data = (on.data * np.exp(-1j * thermal_phase)).real / np.median(
    np.abs(off.data)
)
enhancement.name("enhancement")
# }}}

# {{{ three-line pseudo-Voigt fit (shared center + 14N hyperfine spacing)
nu_offset, E_0, sigma, lambda_L, eta, nu_c, A_hf, A_0, A_1, A_2 = sp.symbols(
    "nu_offset E_0 sigma lambda_L eta nu_c A_hf A_0 A_1 A_2", real=True
)
fitdata = psd.lmfitdata(enhancement)
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
depth = float(abs(enhancement.data).max())
sweep = float(np.ptp(nu))
fitdata.set_guess(
    E_0=dict(value=0.0, min=-depth, max=depth),
    nu_c=dict(value=0.0, min=-10.0, max=10.0),
    A_hf=dict(value=45.0, min=30.0, max=60.0),  # 14N hyperfine spacing
    A_0=dict(value=-depth, min=-20 * depth, max=20 * depth),
    A_1=dict(value=-depth, min=-20 * depth, max=20 * depth),
    A_2=dict(value=-depth, min=-20 * depth, max=20 * depth),
    sigma=dict(value=sweep / 12, min=0.5, max=sweep / 2),
    lambda_L=dict(value=sweep / 12, min=0.5, max=sweep / 2),
    eta=dict(value=0.5, min=0.01, max=0.99),
)
# stage 1: line positions frozen; stage 2: released from the stage-1 result
for guess_name, guess_par in fitdata.guess_parameters.items():
    guess_par.vary = guess_name not in ("nu_c", "A_hf")
fitdata.fit(use_jacobian=False)
fitdata.guess_parameters = fitdata.fit_parameters
for guess_par in fitdata.guess_parameters.values():
    guess_par.vary = True
fitdata.fit(use_jacobian=False)
# }}}

with psd.figlist_var() as fl:
    fl.basename = thisfile
    fl.next("hydroxyTEMPO ODNP field sweep")
    psd.plot(enhancement, "o", color="k", markersize=3, label="DNP curve")
    psd.plot(
        fitdata.eval(500), color="r", alpha=0.8, label="3-line pseudo-Voigt"
    )
    nu_c_fit, A_hf_fit = fitdata.output("nu_c"), fitdata.output("A_hf")
    for this_center in (nu_c_fit - A_hf_fit, nu_c_fit, nu_c_fit + A_hf_fit):
        plt.axvline(this_center, ls=":", color="grey", alpha=0.5)
    plt.axhline(0, lw=0.6, color="k")
    plt.xlabel(r"EPR offset $\nu-\nu_\mathrm{center}$ / MHz")
    plt.ylabel(r"enhancement $S_\mathrm{on}/S_\mathrm{off}$")
    plt.legend()
    plt.gca().grid(alpha=0.3)

    print(f"center field: {center_field_G:#0.8g} G")
    print(f"deepest enhancement: {enhancement.data.min():#0.4g}")
    print(
        f"14N hyperfine spacing: {A_hf_fit:#0.4g} MHz;  centers (MHz): "
        f"{nu_c_fit - A_hf_fit:#0.4g}, {nu_c_fit:#0.4g}, "
        f"{nu_c_fit + A_hf_fit:#0.4g}"
    )
