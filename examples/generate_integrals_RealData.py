"""Generate ODNP integral tables from real data.

This is the updated version of
``examples/broken/generate_integrals_RealData.py``.  It uses the current
``table_of_integrals`` workflow and writes the top-level HDF5 nodes that
``fit_ODNP_data.py`` expects in ``<source filename>_integrals.h5``:

* ``Ep``: normalized enhancement integrals vs microwave power
* ``R1p``: fitted relaxation rates vs microwave power, after FIR correlation
  alignment
* ``T1p``: reciprocal of ``R1p``

Each saved node also carries ``source_file`` and ``Apodization`` metadata so
the downstream fit can identify how the table was generated.
"""

from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd
import sympy as sp
from pyspecProcScripts.generate_coordinates_from_log import (
    generate_coordinates_from_log,
)


plt.rcParams["image.aspect"] = "auto"

# {{{ user-editable block
thisfile, thisexptype, nodename = (
    "260625_hydroxytempo_ODNP_5.h5",
    "B27/ODNP",
    "ODNP",
)
output_dir = Path("/Users/atahan/exp_data/Atahan_Processed_Data/ODNP")
dataset_id = thisfile.removesuffix(".h5")
output_file = f"{dataset_id}_integrals.h5"
output_path = output_dir / output_file
show_alignment_diagnostics = False
alignment_mask_sigma = 150.0
alignment_max_shift = 2000.0
Ep_alignment_max_shift = 1500.0
equal_energy_apodization = True
# }}}


# {{{ Find FIR HDF5 nodes and convert their names to microwave powers
# The real-data files store each inversion-recovery power point as a separate
# top-level node named like ``FIR_12dBm`` or ``FIR_noPower``.  We read those
# node names once from the HDF5 file, parse the dBm value from each name, convert
# it to watts, and sort the resulting list by power so every table written below
# uses a consistent power axis.
filename = psd.search_filename(thisfile, exp_type=thisexptype, unique=True)
with h5py.File(filename, "r") as h5file:
    fir_node_names = [k for k in h5file.keys() if k.startswith("FIR_")]
fir_nodes = []
for thisnodename in fir_node_names:
    m = re.search(r"([0-9]+(?:[p.][0-9]+)?)dBm", thisnodename)
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
    # The enhancement data is loaded exactly like the usual ODNP processing
    # script: attach log data when needed, rebuild the structured indirect
    # coordinate for old post-processing versions, integrate the signal, then
    # normalize all enhancement integrals by the no-power point.  The indirect
    # axis starts as a structured coordinate with both time and power fields;
    # table_of_integrals needs the time coordinate, while the saved table
    # needs the power coordinate.
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
        max_shift=Ep_alignment_max_shift,
        center_aligned_peak=False,
        show_alignment_diagnostics=show_alignment_diagnostics,
        require_unitary_error=False,
        equal_energy_apodization=equal_energy_apodization,
    )
    Ep /= Ep["indirect", 0:1]
    Ep["indirect"] = orig_axis
    Ep.set_error("indirect", orig_axis_error)
    Ep.set_error("indirect", Ep.get_error("indirect")["power"])
    Ep["indirect"] = Ep["indirect"]["power"]
    Ep.set_units("indirect", "W").rename("indirect", "power")
    acq_params = Ep.get_prop("acq_params")
    Ep_to_save = Ep.C
    Ep_to_save.name("Ep")
    Ep_to_save.set_prop("acq_params", acq_params)
    Ep_to_save.set_prop("source_file", thisfile)
    Ep_to_save.set_prop("Apodization", bool(equal_energy_apodization))
    Ep_to_save.hdf5_write(
        output_file,
        directory=str(output_dir),
    )
    print(f"saved Ep -> {output_dir / output_file}")
    # }}}

    # {{{ Prepare node-by-node FIR range detection
    # Each FIR power can have a slightly different apparent line shape, so do
    # not force every node to use the range detected from the highest-power
    # experiment.  Inside the power-ordered loop below, find_peakrange is tried
    # on each raw FIR node before integration.  If a node has multiple apparent
    # peaks or otherwise fails the single-peak test, reuse the most recent
    # successfully detected range from the previous power point.
    previous_signal_range = None
    for thisnodename, _ in fir_nodes:
        s_seed = psd.find_file(
            thisfile,
            exp_type=thisexptype,
            expno=thisnodename,
            lookup=prscr.lookup_table,
        )
        if "nScans" in s_seed.dimlabels:
            s_seed = prscr.clock_correct(s_seed)
        s_seed = s_seed.squeeze()
        try:
            signal_pathway = s_seed.get_prop("coherence_pathway")
            pathway_data = prscr.select_pathway(s_seed.C, signal_pathway)
            frq_center, frq_half = prscr.find_peakrange(
                pathway_data,
                direct="t2",
                peak_lower_thresh=0.05,
            )
            previous_signal_range = tuple(
                sorted((frq_center - frq_half, frq_center + frq_half))
            )
            print(
                "Seeded initial FIR fallback range from"
                f" {thisnodename}: {previous_signal_range}"
            )
            break
        except ValueError:
            pass
    if previous_signal_range is None:
        raise ValueError("could not determine a single FIR signal range")
    # }}}

    # {{{ Fit aligned FIR nodes to generate R1(p), then save T1(p)=1/R1(p)
    # Every FIR node is loaded in power order and aligned before integration.
    # table_of_integrals handles the correlation alignment and integration in
    # one pass.  The integrated inversion-recovery points are fit to the
    # finite-repetition
    # saturation-recovery expression used in the old example.  R1(p) is saved
    # after every power point so a long run leaves a usable partial table on
    # disk, and T1(p) is saved as the reciprocal table expected by the fitting
    # script.
    R1p = psd.ndshape([("power", len(fir_nodes))]).alloc(dtype=np.float64)
    R1p.data[:] = np.nan
    R1p_error = np.nan * np.ones(len(fir_nodes))
    R1p_upper_bound = np.nan * np.ones(len(fir_nodes))
    R1p.set_error(R1p_error.copy())
    R1p.setaxis("power", [pw for _, pw in fir_nodes]).set_units("power", "W")
    R1p.set_prop("acq_params", acq_params)

    for j, (thisnodename, _) in enumerate(fir_nodes):
        fl.basename = thisnodename
        s = psd.find_file(
            thisfile,
            exp_type=thisexptype,
            expno=thisnodename,
            lookup=prscr.lookup_table,
        )
        if "nScans" in s.dimlabels:
            s = prscr.clock_correct(s)
        s = s.squeeze()
        # {{{ Determine this FIR node's frequency range before integration
        # The peak range is still determined in the frequency-domain selected
        # pathway, matching the rough_table_of_integrals style.  The new
        # table_of_integrals function will apply the time-domain apodization,
        # phase correction, correlation alignment, and final integration after
        # receiving this node-specific range.
        try:
            signal_pathway = s.get_prop("coherence_pathway")
            pathway_data = prscr.select_pathway(s.C, signal_pathway)
            frq_center, frq_half = prscr.find_peakrange(
                pathway_data,
                direct="t2",
                peak_lower_thresh=0.05,
            )
            signal_range = tuple(
                sorted((frq_center - frq_half, frq_center + frq_half))
            )
            previous_signal_range = signal_range
            print(
                f"Using FIR signal range from {thisnodename}: {signal_range}"
            )
        except ValueError:
            signal_range = previous_signal_range
            print(
                f"Using previous FIR signal range for {thisnodename}:"
                f" {signal_range}"
            )
        # }}}
        s, ax_last = prscr.table_of_integrals(
            s,
            fl=fl,
            signal_range=signal_range,
            repeat_dims="vd",
            alignment_sigma=alignment_mask_sigma,
            max_shift=alignment_max_shift,
            show_alignment_diagnostics=show_alignment_diagnostics,
            require_unitary_error=False,
            equal_energy_apodization=equal_energy_apodization,
        )
        # {{{ Fit this aligned FIR integral trace to get one R1 value
        # The fit model uses the inversion-recovery recovery delay axis and the
        # finite repetition delay W from acquisition parameters.  This mirrors
        # the old script's methodology while using the current aligned data.
        Mi, R1, vd = sp.symbols("M_inf R_1 vd", real=True)
        acq = s.get_prop("acq_params")
        FIR_rep = acq["FIR_rep"] if "FIR_rep" in acq else acq["FIR_rep_us"]
        W = FIR_rep * 1e-6 + acq["acq_time_ms"] * 1e-3
        vd_units = s.get_units("vd")
        prefactor_scaling = (
            10 ** psd.det_unit_prefactor(vd_units)
            if vd_units is not None
            else 1.0
        )
        R1p_upper_bound[j] = 100 * prefactor_scaling
        f = psd.lmfitdata(s)
        f.functional_form = Mi * (1 - (2 - sp.exp(-W * R1)) * sp.exp(-vd * R1))
        f.set_guess(
            M_inf=dict(
                value=s.max().item(),
                min=0.1 * s.max().item(),
                max=1.5 * s.max().item(),
            ),
            R_1=dict(
                value=0.8 * prefactor_scaling,
                min=0.01 * prefactor_scaling,
                max=R1p_upper_bound[j],
            ),
        )
        f.fit()
        fit_curve = f.eval(200)
        psd.plot(fit_curve, ax=ax_last, alpha=0.5)
        ax_last.text(
            0.5,
            0.5,
            f"RESULT: {f.latex()}",
            ha="center",
            va="center",
            color=fit_curve.get_plot_color(),
            transform=ax_last.transAxes,
        )
        R1p.data[j] = f.output("R_1")
        # {{{ Propagate FIR integral errors through the R1 fit
        # table_of_integrals attaches an uncertainty to each integrated FIR
        # point. I added this because analytical_covariance function in
        # pyspecdata has not been updated. Linearize the fitted inversion-
        # recovery model at the best-fit parameters, build the weighted normal
        # matrix J^T Sigma^-1 J, and take the R1 diagonal element of its
        # inverse as the propagated variance.
        integral_error = s.get_error()
        R1p_error[j] = np.nan
        if integral_error is not None:
            vd_axis = s.getaxis("vd")
            sigma = np.asarray(integral_error, dtype=float)
            finite = np.isfinite(s.data) & np.isfinite(sigma) & (sigma > 0)
            if finite.sum() >= 2:
                M_inf_val = f.output("M_inf")
                R1_val = f.output("R_1")
                exp_W = np.exp(-W * R1_val)
                exp_vd = np.exp(-vd_axis * R1_val)
                dM_inf = 1 - (2 - exp_W) * exp_vd
                dR1 = M_inf_val * exp_vd * ((2 - exp_W) * vd_axis - W * exp_W)
                jacobian = np.vstack([dM_inf, dR1]).T[finite]
                weights = 1.0 / sigma[finite] ** 2
                normal_matrix = jacobian.T @ (weights[:, None] * jacobian)
                covariance = np.linalg.pinv(normal_matrix)
                R1p_error[j] = np.sqrt(covariance[1, 1])
        if not np.isfinite(R1p_error[j]):
            fit_stderr = f.fit_output.params["R_1"].stderr
            if fit_stderr is not None and np.isfinite(fit_stderr):
                R1p_error[j] = fit_stderr
        R1p.set_error(R1p_error.copy())
        # }}}
        # }}}

        # {{{ Save the partial R1(p) table and its reciprocal T1(p) table
        # These writes intentionally happen inside the loop.  If the script is
        # stopped midway through a detailed power series, the HDF5 output still
        # contains all completed powers and NaNs for powers not fit yet.
        R1p_to_save = R1p.C
        R1p_to_save.name("R1p")
        R1p_to_save.set_prop("acq_params", acq_params)
        R1p_to_save.set_prop("source_file", thisfile)
        R1p_to_save.set_prop("Apodization", bool(equal_energy_apodization))
        R1p_to_save.set_prop(
            "R1_fit_upper_bound", float(np.nanmax(R1p_upper_bound))
        )
        R1p_to_save.hdf5_write(
            output_file,
            directory=str(output_dir),
        )
        print(f"saved R1p -> {output_dir / output_file}")
        T1p_to_save = (1.0 / R1p).C
        if R1p.get_error() is not None:
            T1p_to_save.set_error(R1p.get_error() / R1p.data**2)
        T1p_to_save.name("T1p")
        T1p_to_save.set_prop("acq_params", acq_params)
        T1p_to_save.set_prop("source_file", thisfile)
        T1p_to_save.set_prop("Apodization", bool(equal_energy_apodization))
        T1p_to_save.set_prop(
            "R1_fit_upper_bound", float(np.nanmax(R1p_upper_bound))
        )
        T1p_to_save.hdf5_write(
            output_file,
            directory=str(output_dir),
        )
        print(f"saved T1p -> {output_dir / output_file}")
        # }}}
    # }}}
