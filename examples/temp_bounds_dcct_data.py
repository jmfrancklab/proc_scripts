"""Shared real DCCT data for temporary frequency-bound examples."""

import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

import pyspecProcScripts as pypcs
from pyspecProcScripts.load_data import lookup_table

EXPERIMENT_TYPES = ("FIR", "FIR no power", "ODNP enhancement")


def load_dcct_dataset(experiment_type):
    """Load one experimental DCCT dataset and impose known field shifts.

    Both powered and no-power FIR nodes contain a 2x2 ``ph1``/``ph2`` cycle.
    The ODNP node was
    acquired with one four-step ``ph1`` cycle, which likewise contains one
    signal pathway and three artifact pathways.  The stored phase-cycle
    structure is preserved rather than reshaped into dimensions that were not
    acquired.
    """
    if experiment_type in ("FIR", "FIR no power"):
        node_name = (
            "FIR_19.6dBm" if experiment_type == "FIR" else "FIR_noPower"
        )
        postproc = "spincore_IR_v4"
        indirect = "vd"
        seed = 2
    elif experiment_type == "ODNP enhancement":
        node_name = "ODNP"
        postproc = "spincore_ODNP_v6"
        indirect = "indirect"
        seed = 37
    else:
        raise ValueError(f"unknown experiment type {experiment_type!r}")

    raw_data = psd.find_file(
        "260625_hydroxytempo_ODNP_5.h5",
        exp_type="B27/ODNP",
        expno=node_name,
        postproc=postproc,
        lookup=lookup_table,
        print_result=False,
    )
    if experiment_type == "ODNP enhancement":
        raw_data.setaxis(
            "indirect",
            raw_data.getaxis("indirect")["power"],
        ).set_units("indirect", "W").rename("indirect", "power")
        indirect = "power"
    phase_dimensions = [
        dimension
        for dimension in raw_data.dimlabels
        if dimension.startswith("ph")
    ]
    max_alignment_shift = 450.0
    shift_dimensions = phase_dimensions + [indirect]
    shift_shape = tuple(
        raw_data.shape[dimension] for dimension in shift_dimensions
    )
    rng = np.random.default_rng(seed)
    frequency_shifts = psd.nddata(
        np.clip(
            rng.uniform(
                -0.85 * max_alignment_shift,
                0.85 * max_alignment_shift,
                size=(1,) * len(phase_dimensions)
                + (raw_data.shape[indirect],),
            )
            + rng.normal(scale=20.0, size=shift_shape),
            -0.9 * max_alignment_shift,
            0.9 * max_alignment_shift,
        ),
        shift_shape,
        shift_dimensions,
    )
    for dimension in shift_dimensions:
        frequency_shifts.setaxis(dimension, raw_data.getaxis(dimension))

    # Isolate the measured signal pathway so that the other three pathways in
    # the controlled comparison arise specifically from the imposed field
    # instability rather than artifacts already present in the acquisition.
    shifted_data = raw_data.C
    for dimension, coherence in raw_data.get_prop("coherence_pathway").items():
        shifted_data *= shifted_data.fromaxis(dimension).run(
            lambda coordinate, selected=coherence: coordinate == selected
        )

    # A field offset acts on individual phase-cycle transients.  Transforming
    # the shifted data back to the coherence domain transfers some of the
    # signal into the three nominally empty pathways, as in a real unstable
    # field acquisition.
    shifted_data.ift("t2").ift(phase_dimensions)
    shifted_data *= np.exp(
        1j * 2 * np.pi * frequency_shifts * shifted_data.fromaxis("t2")
    )
    shifted_data.ft("t2").ft(phase_dimensions)
    return (
        raw_data,
        shifted_data,
        {
            "experiment_type": experiment_type,
            "indirect": indirect,
            "frequency_shifts": frequency_shifts,
            "phase_dimensions": phase_dimensions,
            "signal_pathway": raw_data.get_prop("coherence_pathway"),
            "max_alignment_shift": max_alignment_shift,
        },
    )


def show_raw_data(fl, raw_data, configuration):
    """Show the untouched experimental DCCT data."""
    fl.next(f"{configuration['experiment_type']}: raw experimental data")
    psd.DCCT(
        raw_data,
        fig=plt.gcf(),
        title="raw experimental data",
    )


def sign_corrected(data, configuration, inhomogeneous_bounds):
    """Return a copy with FIR or enhancement sign changes removed."""
    signal = pypcs.select_pathway(
        data,
        configuration["signal_pathway"],
    )
    signs = pypcs.determine_sign(signal, inhomogeneous_bounds)
    corrected = data * signs
    # The demonstrations use no uncertainty-weighted operations, while their
    # zero-filled transforms change the direct-axis length.  Drop the stored
    # pointwise errors so they are not incorrectly reshaped after padding.
    corrected.set_error(None)
    return corrected


def determine_preliminary_bounds(data, configuration):
    """Return echo timing and narrow bounds without expanding the bandwidth.

    Linewidth and centering diagnostics need these bounds for sign correction,
    but do not require an alignment window to fit inside the acquired spectrum.
    Work on a copy so preliminary detection leaves the source unchanged.
    """
    working_data = data.C
    echo_center = pypcs.find_exponential_echo_center(
        working_data,
        decay_rate=250,
    )
    working_data.set_prop("echo_center", echo_center)
    pypcs.det_inh_bounds(
        working_data,
        0.10,
        peak_lowest_thresh=0.03,
        signal_pathway=configuration["signal_pathway"],
    )
    return {
        "echo_center": echo_center,
        "inh_bounds": working_data.get_prop("inh_bounds").copy(),
    }


def determine_demo_bounds(data, configuration):
    """Determine narrow bounds and validate the expanded processing window."""
    bounds = determine_preliminary_bounds(data, configuration)
    working_data = data.C
    working_data.set_prop("echo_center", bounds["echo_center"])
    processed = pypcs.fid_from_echo(
        working_data,
        configuration["signal_pathway"],
        inh_bounds=bounds["inh_bounds"],
        max_alignment_shift=configuration["max_alignment_shift"],
    )
    return {
        **bounds,
        "homogeneous_linewidth": processed.get_prop("homogeneous_linewidth"),
        "processing_bounds": processed.get_prop("processing_bounds").copy(),
    }


# SINGLE_USE_EXCEPTION -- callback factory required by alignment API
def frequency_mask(bounds):
    """Return the square-root frequency mask required by ``correl_align``."""

    def apply_mask(data):
        return data * data.fromaxis("t2").run(
            lambda frequency: (frequency >= bounds[0])
            & (frequency <= bounds[1])
        )

    return apply_mask


# SINGLE_USE_EXCEPTION -- callback required by correlation-alignment API
def keep_all_coherence_pathways(coherence_mask):
    """Use the signal pathway and all three field-shift artifacts."""
    coherence_mask.data[:] = 1
    return coherence_mask
