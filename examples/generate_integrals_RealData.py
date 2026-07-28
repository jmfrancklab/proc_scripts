"""Generate ODNP integral tables from real data.

This script writes the top-level HDF5 nodes expected by the ODNP fitting
script:

* ``Ep``: normalized enhancement integrals vs microwave power
* ``R1p``: fitted relaxation rates vs microwave power
* ``T1p``: reciprocal relaxation times vs microwave power

The FIR integrations use ``table_of_integrals``.  Node discovery and cached
fallback peak ranges live here because they depend on the source HDF5 file
layout, while DC, clock, phase, alignment, and integration processing remain
inside the shared pipeline.
"""

from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
from sympy import exp as s_exp
from sympy import symbols
from pyspecProcScripts.generate_coordinates_from_log import (
    generate_coordinates_from_log,
)

plt.rcParams["image.aspect"] = "auto"
if not hasattr(psd.lmfitdata, "settoguess"):
    psd.lmfitdata.settoguess = psd.lmfitdata.set_to_guess

# {{{ changeable parameters
thisfile, thisexptype, nodename = (
    #    "260724_TMTPDI_ODNP_1.h5",
    "260625_hydroxytempo_ODNP_5.h5",
    "B27/ODNP",
    "ODNP",
)
output_dir = Path("/Users/atahan/exp_data/Atahan_Processed_Data/ODNP")
dataset_id = thisfile.removesuffix(".h5")
output_file = f"{dataset_id}_integrals.h5"
show_alignment_diagnostics = False
alignment_mask_sigma = 125.0
Ep_alignment_max_shift_Hz = 440 * 2.5
fid_from_echo_slice_multiplier = 5
# }}}


# {{{ Find FIR HDF5 nodes and convert their names to microwave powers
filename = psd.search_filename(thisfile, exp_type=thisexptype, unique=True)
with h5py.File(filename, "r") as h5file:
    fir_node_names = [
        name
        for name in h5file.keys()
        if re.match(r"^FIR_(?:noPower|-?\d+(?:[p.]\d+)?dBm)$", name)
    ]
fir_nodes = []
for thisnodename in fir_node_names:
    m = re.search(r"(-?\d+(?:[p.]\d+)?)dBm", thisnodename)
    if m:
        power = prscr.dBm2power(float(m.group(1).replace("p", "."))).item()
    elif "noPower" in thisnodename:
        power = 0.0
    else:
        raise ValueError(f"Cannot infer power from node name {thisnodename!r}")
    fir_nodes.append((thisnodename, power))
fir_nodes.sort(key=lambda x: x[1])
print(f"Found {len(fir_nodes)} FIR node(s): {[n for n, _ in fir_nodes]}")
# }}}

with psd.figlist_var() as fl:
    fl.basename = thisfile

    # {{{ Generate normalized E(p) table from the enhancement node
    Ep = psd.find_file(
        thisfile,
        exp_type=thisexptype,
        expno=nodename,
        lookup=prscr.lookup_table,
        fl=fl,
    )
    if Ep.get_prop("log") is None:
        Ep = prscr.attach_log_data_from_file(Ep, thisfile, thisexptype)
    m = re.search(r".*ODNP.*v([0-9]+)$", Ep.get_prop("postproc_type"))
    if m is None:
        raise IOError(
            f"Unexpected postproc_type: {Ep.get_prop('postproc_type')!r}"
        )
    if int(m.groups()[0]) < 6:
        Ep = generate_coordinates_from_log(Ep, fl=fl)
    orig_axis = Ep["indirect"]
    orig_axis_error = Ep.get_error("indirect")
    Ep["indirect"] = Ep["indirect"]["time"]
    Ep.set_units("indirect", "s")
    Ep.set_error("indirect", orig_axis_error["time"])
    Ep, _ = prscr.table_of_integrals(
        Ep,
        fl=fl,
        repeat_dims="indirect",
        alignment_sigma=alignment_mask_sigma,
        max_shift_Hz=Ep_alignment_max_shift_Hz,
        center_aligned_peak=False,
        show_alignment_diagnostics=show_alignment_diagnostics,
        fid_from_echo_slice_multiplier=fid_from_echo_slice_multiplier,
    )
    Ep["indirect"] = orig_axis
    Ep.set_error("indirect", orig_axis_error)
    Ep.set_error("indirect", Ep.get_error("indirect")["power"])
    Ep["indirect"] = Ep["indirect"]["power"]
    Ep.set_units("indirect", "W").rename("indirect", "power")
    # normalize() propagates covariance for real data, so first phase the
    # complex reference onto the real axis without changing its magnitude.
    Ep_reference = Ep["power", 0].item()
    Ep /= Ep_reference / abs(Ep_reference)
    Ep = Ep.real.normalize("power")
    acq_params = Ep.get_prop("acq_params")
    Ep_to_save = Ep.C
    Ep_to_save.name("Ep")
    Ep_to_save.set_prop("acq_params", acq_params)
    Ep_to_save.set_prop("source_file", thisfile)
    Ep_to_save.hdf5_write(output_file, directory=str(output_dir))
    print(f"saved Ep -> {output_dir / output_file}")
    # }}}

    # {{{ Fit aligned FIR nodes to generate R1(p) and T1(p)
    R1p = psd.ndshape([("power", len(fir_nodes))]).alloc(dtype=np.float64)
    R1p.data[:] = np.nan
    R1p_error = np.nan * np.ones(len(fir_nodes))
    R1p.set_error(R1p_error.copy())
    R1p.setaxis("power", [pw for _, pw in fir_nodes]).set_units("power", "W")
    R1p.set_prop("acq_params", acq_params)

    previous_signal_range = None
    previous_fallback_node = None
    for j, (thisnodename, _) in sorted(
        enumerate(fir_nodes),
        key=lambda indexed_node: indexed_node[1][1],
        reverse=True,
    ):
        fl.basename = thisnodename
        s = psd.find_file(
            thisfile,
            exp_type=thisexptype,
            expno=thisnodename,
            lookup=prscr.lookup_table,
        )
        s = s.squeeze()
        signal_pathway = s.get_prop("coherence_pathway")
        acq = s.get_prop("acq_params")
        align_max_shift_hz = 2.5 * acq["tolerance_Hz"]

        # {{{ Cache this node's raw range for the next lower-power node
        # The current node first searches its own fully processed signal inside
        # table_of_integrals.  Only if that search fails does it use the raw
        # range cached from the preceding higher-power node.
        fallback_signal_range = previous_signal_range
        fallback_node = previous_fallback_node
        current_signal_range = None
        pathway_data = prscr.select_pathway(s.C, signal_pathway)
        pathway_data = pathway_data.C
        for dimname in ("vd", "nScans", "repeats"):
            if dimname in pathway_data.dimlabels:
                pathway_data = pathway_data.mean(dimname)
        try:
            frq_center, frq_half = prscr.find_peakrange(
                pathway_data,
                direct="t2",
                peak_lower_thresh=0.1,
            )
            frq_half = abs(frq_half)
            current_signal_range = tuple(
                sorted(frq_center + np.r_[-1, 1] * frq_half)
            )
            print(
                f"Cached FIR fallback range from {thisnodename}: "
                f"{current_signal_range}"
            )
        except ValueError as e:
            print(f"Could not cache a raw range from {thisnodename} ({e})")
        # }}}

        s, ax_last = prscr.table_of_integrals(
            s,
            fl=fl,
            signal_pathway=signal_pathway,
            repeat_dims="vd",
            alignment_sigma=alignment_mask_sigma,
            max_shift_Hz=align_max_shift_hz,
            center_aligned_peak=True,
            show_alignment_diagnostics=show_alignment_diagnostics,
            fallback_signal_range=fallback_signal_range,
            clock_correction=True,
            fid_from_echo_slice_multiplier=fid_from_echo_slice_multiplier,
        )
        if s.get_prop("table_of_integrals_used_fallback"):
            print(
                f"{thisnodename} used the remembered range from "
                f"{fallback_node}: {fallback_signal_range}"
            )
        if current_signal_range is not None:
            previous_signal_range = current_signal_range
            previous_fallback_node = thisnodename

        # {{{ Fit this aligned FIR integral trace to get one R1 value
        fit_data = s.C.run(np.real)
        M_inf, R_1, vd = symbols("M_inf R_1 vd", real=True)
        if "FIR_rep" in acq:
            repetition_s = acq["FIR_rep"] * 1e-6
        elif "FIR_rep_us" in acq:
            repetition_s = acq["FIR_rep_us"] * 1e-6
        else:
            raise KeyError("acq_params must contain FIR_rep or FIR_rep_us")
        W = repetition_s + acq.get("acq_time_ms", 0) * 1e-3
        vd_units = s.get_units("vd")
        prefactor_scaling = (
            10 ** psd.det_unit_prefactor(vd_units)
            if vd_units is not None
            else 1.0
        )
        r1_estimates_per_s = []
        if repetition_s > 0:
            r1_estimates_per_s.append(2.0 / repetition_s)
        if all(k in acq for k in ("concentration", "krho_hot", "T1water_hot")):
            r1_estimates_per_s.append(
                acq["concentration"] * acq["krho_hot"]
                + 1.0 / acq["T1water_hot"]
            )
        if all(
            k in acq for k in ("concentration", "krho_cold", "T1water_cold")
        ):
            r1_estimates_per_s.append(
                acq["concentration"] * acq["krho_cold"]
                + 1.0 / acq["T1water_cold"]
            )
        finite_r1_estimates = [
            float(x) for x in r1_estimates_per_s if np.isfinite(x) and x > 0
        ]
        if finite_r1_estimates:
            r1_guess_per_s = float(np.median(finite_r1_estimates))
        else:
            positive_vd = np.asarray(s.getaxis("vd"), dtype=float)
            positive_vd = positive_vd[
                np.isfinite(positive_vd) & (positive_vd > 0)
            ]
            if positive_vd.size == 0:
                raise ValueError(
                    "Cannot infer an R1 guess without positive vd values"
                )
            r1_guess_per_s = 1.0 / np.median(positive_vd)
        r1_bounds_per_s = (
            max(1e-6, r1_guess_per_s / 5.0),
            r1_guess_per_s * 5.0,
        )
        y_data = np.asarray(fit_data.data, dtype=float)
        finite_y = y_data[np.isfinite(y_data)]
        if finite_y.size == 0:
            raise ValueError(
                "Cannot infer M_inf bounds from non-finite FIR integrals"
            )
        signal_scale = max(float(np.nanmax(np.abs(finite_y))), 1.0)
        signal_guess = float(
            np.nanmedian(finite_y[-max(1, finite_y.size // 3) :])
        )
        if abs(signal_guess) < 0.05 * signal_scale:
            signal_guess = float(finite_y[np.nanargmax(np.abs(finite_y))])
        if signal_guess >= 0:
            M_inf_guess = max(signal_guess, 0.1 * signal_scale)
            M_inf_bounds = (0.0, 2.0 * signal_scale)
        else:
            M_inf_guess = min(signal_guess, -0.1 * signal_scale)
            M_inf_bounds = (-2.0 * signal_scale, 0.0)
        fit = psd.lmfitdata(fit_data)
        fit.functional_form = M_inf * (
            1 - (2 - s_exp(-W * R_1)) * s_exp(-vd * R_1)
        )
        fit.set_guess(
            M_inf=dict(
                value=M_inf_guess,
                min=M_inf_bounds[0],
                max=M_inf_bounds[1],
            ),
            R_1=dict(
                value=r1_guess_per_s * prefactor_scaling,
                min=r1_bounds_per_s[0] * prefactor_scaling,
                max=r1_bounds_per_s[1] * prefactor_scaling,
            ),
        )
        fit.settoguess()
        guess = fit.eval(100)
        fit.fit()
        fit_curve = fit.eval(100)
        R1 = fit.output("R_1") / prefactor_scaling
        R1_stderr = fit.fit_output.params["R_1"].stderr
        R1_uncertainty = (
            R1_stderr / prefactor_scaling
            if R1_stderr is not None and np.isfinite(R1_stderr)
            else np.nan
        )
        R1p.data[j] = R1
        R1p_error[j] = R1_uncertainty
        R1p.set_error(R1p_error.copy())
        psd.plot(guess, "-", ax=ax_last, alpha=0.25, human_units=False)
        psd.plot(fit_curve, "-", ax=ax_last, alpha=0.7, human_units=False)
        ax_last.text(
            0.5,
            0.5,
            (
                rf"$R_1={R1:.3g}\pm{R1_uncertainty:.2g}"
                r"\ \mathrm{s^{-1}}$"
                f"\nT1={1.0 / R1:.3g} s"
            ),
            ha="center",
            va="center",
            color=fit_curve.get_plot_color(),
            transform=ax_last.transAxes,
        )
        print(f"{thisnodename}: R1={R1:#0.6g} s^-1, T1={1.0 / R1:#0.6g} s")
        # }}}

        # {{{ Save partial R1(p) and T1(p) tables after each FIR node
        R1p_to_save = R1p.C
        R1p_to_save.name("R1p")
        R1p_to_save.set_prop("acq_params", acq_params)
        R1p_to_save.set_prop("source_file", thisfile)
        R1p_to_save.hdf5_write(output_file, directory=str(output_dir))
        print(f"saved R1p -> {output_dir / output_file}")
        T1p_to_save = (1.0 / R1p).C
        if R1p.get_error() is not None:
            T1p_to_save.set_error(R1p.get_error() / R1p.data**2)
        T1p_to_save.name("T1p")
        T1p_to_save.set_prop("acq_params", acq_params)
        T1p_to_save.set_prop("source_file", thisfile)
        T1p_to_save.hdf5_write(output_file, directory=str(output_dir))
        print(f"saved T1p -> {output_dir / output_file}")
        # }}}
    # }}}
