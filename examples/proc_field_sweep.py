"""
Process field-swept enhancement data
====================================

``py proc_field_sweep.py NODENAME FILENAME EXP_TYPE``

Process a field sweep acquired by ``run_field_sweep.py``.  The acquisition
stores the requested and measured fields in the ``field_axis_G`` and
``field_readback_G`` properties.  Here we use the measured field readback as
the sweep axis, integrate each echo with ``rough_table_of_integrals``,
normalize to the off-resonance edge points, and fit the field profile with a
simple pseudo-Voigt line shape using ``lmfitdata``.

Tested with:

``py proc_field_sweep.py field_sweep 260000_field_sweep.h5\
        ODNP_NMR_comp/field_dependent``
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp

plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2
plt.rcParams.update(
    {
        "errorbar.capsize": 2,
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.figsize": (6, 4),
    }
)

if (
    "SPHINX_GALLERY_RUNNING" in os.environ
    and os.environ["SPHINX_GALLERY_RUNNING"] == "True"
):
    sys.argv = [
        sys.argv[0],
        "field_sweep",
        "260000_field_sweep.h5",  # TODO: Change the file name once run
        "B27/field_dependent",
    ]

assert len(sys.argv) == 4, "intended to be called with file info at cmdline"
s = psd.find_file(
    sys.argv[2],
    exp_type=sys.argv[3],
    expno=sys.argv[1],
    lookup=prscr.lookup_table,
)

with psd.figlist_var() as fl:
    fl.basename = sys.argv[2]
    acq_params = s.get_prop("acq_params")
    field_G = np.asarray(s.get_prop("field_readback_G"), dtype=float)
    center_field_G = (
        acq_params["carrierFreq_MHz"] / acq_params["gamma_eff_MHz_G"]
    )
    s["indirect"] = (field_G - center_field_G) * acq_params["gamma_eff_MHz_G"]
    s.rename("indirect", "nu_offset")
    s.set_units("nu_offset", "MHz")
    s.sort("nu_offset")
    s, _ = prscr.rough_table_of_integrals(s, fl=fl, title=sys.argv[2])
    if s.get_units("nu_offset") == "kHz":
        s["nu_offset"] = s["nu_offset"] / 1e3
        s.set_units("nu_offset", "MHz")

    # Treat the sweep edges as off-resonance signal, so E(p)=1 there.
    edge_idx = np.r_[0:1, s.shape["nu_offset"] - 1 : s.shape["nu_offset"]]
    edge_signal = s.data.real[edge_idx].mean()
    if np.isclose(edge_signal, 0):
        raise ValueError("Cannot normalize field sweep: edge signal is zero")
    s /= edge_signal
    s.name("enhancement")

    nu_offset, E0, A, nu_center, sigma, lambda_L, eta = sp.symbols(
        "nu_offset E_0 A nu_center sigma lambda_L eta", real=True
    )
    fitdata = psd.lmfitdata(s)
    fitdata.functional_form = E0 + A * (
        eta * lambda_L**2 / ((nu_offset - nu_center) ** 2 + lambda_L**2)
        + (1 - eta) * sp.exp(-((nu_offset - nu_center) ** 2) / (2 * sigma**2))
    )
    baseline_guess = s.data.real[edge_idx].mean()
    y_from_baseline = s.data.real - baseline_guess
    nu_axis = np.asarray(s["nu_offset"], dtype=float)
    sweep_width = nu_axis.max() - nu_axis.min()
    linewidth_guess = max(sweep_width / 8, 1e-3)
    amp_guess = y_from_baseline[np.argmax(abs(y_from_baseline))]
    fitdata.set_guess(
        E_0=dict(value=baseline_guess, min=0.2, max=2.0),
        A=dict(value=amp_guess, min=2 * amp_guess, max=0)
        if amp_guess < 0
        else dict(value=amp_guess, min=0, max=2 * amp_guess),
        nu_center=dict(
            value=nu_axis[np.argmax(abs(y_from_baseline))],
            min=nu_axis.min(),
            max=nu_axis.max(),
        ),
        sigma=dict(
            value=linewidth_guess,
            min=1e-4,
            max=max(sweep_width, 2 * linewidth_guess),
        ),
        lambda_L=dict(
            value=linewidth_guess,
            min=1e-4,
            max=max(sweep_width, 2 * linewidth_guess),
        ),
        eta=dict(value=0.5, min=0, max=1),
    )
    for name, par in fitdata.guess_parameters.items():
        if name != "nu_center":
            par.vary = False
    fitdata.fit(use_jacobian=False)
    fitdata.guess_parameters = fitdata.fit_parameters
    for name, par in fitdata.guess_parameters.items():
        par.vary = True
    fitdata.fit(use_jacobian=False)
    fit = fitdata.eval(400)
    peak = (
        fit.C.run(lambda x: abs(x - fitdata.output("E_0")))
        .argmax("nu_offset")
        .item()
    )
    peak_enhancement = fit["nu_offset":peak].item().real

    fl.next("field-swept $E(p)$", legend=True)
    ax = plt.gca()
    ax.plot(
        np.asarray(s["nu_offset"], dtype=float),
        s.data.real,
        "o",
        label=r"$E(p)$",
    )
    ax.plot(
        np.asarray(fit["nu_offset"], dtype=float),
        fit.data.real,
        label="pseudo-Voigt fit",
        alpha=0.7,
    )
    ax.axvline(peak, ls=":", color="k", alpha=0.5)
    ax.text(
        peak,
        0.95,
        f" {peak:#0.4g} MHz\n E={peak_enhancement:#0.4g}",
        ha="left",
        va="top",
        color="k",
        transform=mpl.transforms.blended_transform_factory(
            ax.transData, ax.transAxes
        ),
    )
    ax.set_xlabel(r"$\nu-\nu_\mathrm{center}$ / MHz")
    ax.set_ylabel(r"$E(p)$")
    ax.grid()

    print(f"center field: {center_field_G:#0.8g} G")
    print(f"edge-reference signal: {edge_signal:#0.8g}")
    print(fitdata.output())
    print(f"peak enhancement: {peak_enhancement:#0.8g} at {peak:#0.8g} MHz")
