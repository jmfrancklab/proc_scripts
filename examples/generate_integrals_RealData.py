"""Generate ODNP integral tables from real data.

This is the updated version of
``examples/broken/generate_integrals_RealData.py``.  It uses the current
``rough_table_of_integrals`` workflow and writes the tables that
``fit_ODNP_data.py`` expects:

* ``Ep``: normalized enhancement integrals vs microwave power
* ``R1p``: fitted relaxation rates vs microwave power
* ``T1p``: reciprocal of ``R1p``
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
    "260615_hydroxytempo_ODNP_3.h5",
    "B27/ODNP",
    "ODNP",
)
output_dir = Path(__file__).resolve().parent / "generated"
output_file = "ODNP_integrals_27mM.h5"
dataset_id = Path(thisfile).stem
fir_range_expansion = 1.0
# }}}


def node_power(nodename):
    m = re.search(r"([0-9]+(?:p[0-9]+)?)dBm", nodename)
    if m:
        return prscr.dBm2power(float(m.group(1).replace("p", "."))).item()
    if "noPower" in nodename:
        return 0.0
    raise ValueError(f"Cannot infer power from node name {nodename!r}")


def sorted_fir_nodes(h5file, exp_type):
    filename = psd.search_filename(h5file, exp_type=exp_type, unique=True)
    with h5py.File(filename, "r") as f:
        nodes = [k for k in f.keys() if k.startswith("FIR_")]
    return sorted([(n, node_power(n)) for n in nodes], key=lambda x: x[1])


def save_table(table, table_name, acq_params):
    output_dir.mkdir(exist_ok=True)
    table = table.C
    table.name(table_name)
    table.set_prop("acq_params", acq_params)
    table.hdf5_write(f"{output_file}/{dataset_id}", directory=str(output_dir))
    print(f"saved {dataset_id}/{table_name} -> {output_dir / output_file}")


def load_dataset(expno, fl=None):
    return psd.find_file(
        thisfile,
        exp_type=thisexptype,
        expno=expno,
        lookup=prscr.lookup_table,
        fl=fl,
    )


def process_ep(fl):
    s = load_dataset(nodename, fl=fl)
    if s.get_prop("log") is None:
        s = prscr.attach_log_data_from_file(s, thisfile, thisexptype)
    m = re.search(r".*ODNP.*v([0-9]+)$", s.get_prop("postproc_type"))
    if m is None:
        raise IOError(f"Unexpected postproc_type: {s.get_prop('postproc_type')!r}")
    if int(m.groups()[0]) < 6:
        s = generate_coordinates_from_log(s, fl=fl)
    orig_axis = s["indirect"]
    orig_axis_error = s.get_error("indirect")
    s["indirect"] = s["indirect"]["time"]
    s.set_units("indirect", "s")
    s.set_error("indirect", orig_axis_error["time"])
    s, _ = prscr.rough_table_of_integrals(s, fl=fl)
    s.set_error(s["indirect", 0].item() * 0.01)
    s /= s["indirect", 0:1]
    s["indirect"] = orig_axis
    s.set_error("indirect", orig_axis_error)
    s.set_error("indirect", s.get_error("indirect")["power"])
    s["indirect"] = s["indirect"]["power"]
    s.set_units("indirect", "W").rename("indirect", "power")
    return s


def fir_signal_range(s, direct="t2", peak_lower_thresh=0.05):
    signal_pathway = s.get_prop("coherence_pathway")
    pathway_data = prscr.select_pathway(s.C, signal_pathway)
    frq_center, frq_half = prscr.find_peakrange(
        pathway_data, direct=direct, peak_lower_thresh=peak_lower_thresh
    )
    return tuple(sorted((frq_center - frq_half, frq_center + frq_half)))


def expand_range(signal_range, expansion):
    center = np.mean(signal_range)
    half_width = 0.5 * np.diff(signal_range).item() * expansion
    return (center - half_width, center + half_width)


def infer_fir_signal_range(fir_nodes):
    for thisnodename, _ in reversed(fir_nodes):
        s = load_dataset(thisnodename)
        if "nScans" in s.dimlabels:
            s = prscr.clock_correct(s)
        s = s.squeeze()
        try:
            signal_range = expand_range(
                fir_signal_range(s),
                fir_range_expansion,
            )
            print(f"Using FIR signal range from {thisnodename}: {signal_range}")
            return signal_range
        except ValueError:
            pass
    raise ValueError("could not determine a single FIR signal range")


def fit_fir_node(s, ax_last):
    Mi, R1, vd = sp.symbols("M_inf R_1 vd", real=True)
    acq = s.get_prop("acq_params")
    FIR_rep = acq["FIR_rep"] if "FIR_rep" in acq else acq["FIR_rep_us"]
    W = FIR_rep * 1e-6 + acq["acq_time_ms"] * 1e-3
    vd_units = s.get_units("vd")
    prefactor_scaling = (
        10 ** psd.det_unit_prefactor(vd_units) if vd_units is not None else 1.0
    )
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
            max=100 * prefactor_scaling,
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
    return f.output("R_1")


fir_nodes = sorted_fir_nodes(thisfile, thisexptype)
print(f"Found {len(fir_nodes)} FIR node(s): {[n for n, _ in fir_nodes]}")

with psd.figlist_var() as fl:
    fl.basename = thisfile
    Ep = process_ep(fl)
    acq_params = Ep.get_prop("acq_params")
    save_table(Ep, "Ep", acq_params)

    signal_range = infer_fir_signal_range(fir_nodes)
    R1p = psd.ndshape([("power", len(fir_nodes))]).alloc(dtype=np.float64)
    R1p.data[:] = np.nan
    R1p.setaxis("power", [pw for _, pw in fir_nodes]).set_units("power", "W")
    R1p.set_prop("acq_params", acq_params)

    for j, (thisnodename, _) in enumerate(fir_nodes):
        fl.basename = thisnodename
        s = load_dataset(thisnodename)
        if "nScans" in s.dimlabels:
            s = prscr.clock_correct(s)
        s = s.squeeze()
        s, ax_last = prscr.rough_table_of_integrals(
            s,
            fl=fl,
            signal_range=signal_range,
        )
        R1p["power", j] = fit_fir_node(s, ax_last)
        save_table(R1p, "R1p", acq_params)
        save_table(1.0 / R1p, "T1p", acq_params)
