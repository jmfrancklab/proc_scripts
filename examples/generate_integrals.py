"""Show synthetic and real-data integral/alignment examples.

This is the updated version of ``examples/broken/generate_integrals.py``.
The first two examples generate synthetic inversion-recovery and enhancement
data, then compare ``table_of_integrals`` against
``rough_table_of_integrals`` on the same fake data.
The real-data section is diagnostic-only: it compares the correlation
alignment path in ``table_of_integrals`` with the older
``rough_table_of_integrals`` alignment for one FIR node and for E(p).  It does
not write any HDF5 output.
"""

from collections import OrderedDict
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


np.random.seed(2021)
plt.rcParams["image.aspect"] = "auto"

t2, vd, power, ph1, ph2 = sp.symbols("t2 vd power ph1 ph2")
echo_time = 10e-3

examples = [
    dict(
        label="IR",
        expression=23
        * (1 - 2 * sp.exp(-vd / 0.2))
        * sp.exp(1j * 2 * sp.pi * 100 * t2 - abs(t2) * 50 * sp.pi),
        variable_defs=OrderedDict(
            [
                ("vd", psd.nddata(np.r_[0:1:40j], "vd")),
                ("ph1", psd.nddata(np.r_[0:4] / 4.0, "ph1")),
                ("ph2", psd.nddata(np.r_[0, 2] / 4.0, "ph2")),
                ("t2", psd.nddata(np.r_[0:0.2:256j] - echo_time, "t2")),
            ]
        ),
        signal_pathway={"ph1": 0, "ph2": 1},
        indirect="vd",
        signal_range=(-400, 400),
    ),
    dict(
        label="Enhancement",
        expression=23
        * (1 - (32 * power / (0.25 + power)) * 150e-6 * 659.33)
        * sp.exp(1j * 2 * sp.pi * 100 * t2 - abs(t2) * 50 * sp.pi),
        variable_defs=OrderedDict(
            [
                ("power", psd.nddata(np.r_[0:4:25j], "power")),
                ("ph1", psd.nddata(np.r_[0:4] / 4.0, "ph1")),
                ("t2", psd.nddata(np.r_[0:0.2:256j] - echo_time, "t2")),
            ]
        ),
        signal_pathway={"ph1": 1},
        indirect="power",
        signal_range=(-200, 600),
    ),
]


with psd.figlist_var() as fl:
    for cfg in examples:
        fl.basename = f"(fake {cfg['label']})"
        data = psd.fake_data(
            cfg["expression"],
            cfg["variable_defs"],
            cfg["signal_pathway"],
        )
        data.reorder([cfg["indirect"], "t2"], first=False)
        data.ft("t2")
        data /= np.sqrt(psd.ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        corr_int, _ = prscr.table_of_integrals(
            data.C,
            signal_range=cfg["signal_range"],
            signal_pathway=cfg["signal_pathway"],
            fl=fl,
            title=f"fake {cfg['label']} correlation alignment",
            repeat_dims=cfg["indirect"],
            center_aligned_peak=cfg["label"] != "Enhancement",
            equal_energy_apodization=False,
        )
        rough_int, _ = prscr.rough_table_of_integrals(
            data.C,
            signal_range=cfg["signal_range"],
            signal_pathway=cfg["signal_pathway"],
            fl=fl,
            title=f"fake {cfg['label']} rough alignment",
        )
        fl.next(f"fake {cfg['label']} table-of-integrals comparison")
        fl.plot(
            corr_int.C.set_plot_color("k"),
            "o",
            label="fake table_of_integrals correlation alignment",
        )
        fl.plot(
            rough_int.C.set_plot_color("r"),
            "x",
            label="fake rough_table_of_integrals alignment",
        )
        plt.gca().set_title(
            f"fake {cfg['label']} table-of-integrals comparison"
        )
        plt.legend()

    # {{{ Show real-data correlation alignment and rough-table alignment
    # These examples use the most recent 27 mM hydroxytempo ODNP dataset.  The
    # FIR node behaves like a repeated FID/inversion-recovery experiment at one
    # microwave power, so the aligned average is a good representative spectrum
    # and table_of_integrals is allowed to recenter the final window on the
    # aligned peak.  E(p) is different: the indirect dimension is a power
    # trajectory, not repeated measurements of one spectrum.  High-power points
    # can dominate the average, so recentering can move the integration window
    # and introduce sign/range artifacts; therefore center_aligned_peak is
    # disabled for E(p).
    #
    # Nothing is saved here.  We call the two table functions only to show their
    # diagnostic alignment figures side-by-side in the figure list.
    real_file, real_exptype, ep_node = (
        "260622_hydroxytempo_ODNP_4.h5",
        "B27/ODNP",
        "ODNP",
    )
    show_alignment_diagnostics = True
    alignment_mask_sigma = 150.0
    alignment_max_shift = 2000.0
    Ep_alignment_max_shift = 1500.0
    equal_energy_apodization = True
    filename = psd.search_filename(
        real_file, exp_type=real_exptype, unique=True
    )
    fir_node = "FIR_35dBm"

    fl.basename = f"{real_file} {fir_node}"
    fid_data = psd.find_file(
        real_file,
        exp_type=real_exptype,
        expno=fir_node,
        lookup=prscr.lookup_table,
    )
    if "nScans" in fid_data.dimlabels:
        fid_data = prscr.clock_correct(fid_data)
    fid_data = fid_data.squeeze()
    fid_pathway = fid_data.get_prop("coherence_pathway")
    frq_center, frq_half = prscr.find_peakrange(
        prscr.select_pathway(fid_data.C, fid_pathway),
        direct="t2",
        peak_lower_thresh=0.05,
    )
    fid_signal_range = tuple(
        sorted((frq_center - frq_half, frq_center + frq_half))
    )
    prscr.table_of_integrals(
        fid_data.C,
        signal_range=fid_signal_range,
        signal_pathway=fid_pathway,
        fl=fl,
        title=f"{fir_node} correlation alignment",
        repeat_dims="vd",
        alignment_sigma=alignment_mask_sigma,
        max_shift=alignment_max_shift,
        center_aligned_peak=True,
        show_alignment_diagnostics=show_alignment_diagnostics,
        require_unitary_error=False,
        equal_energy_apodization=equal_energy_apodization,
    )
    prscr.rough_table_of_integrals(
        fid_data.C,
        signal_range=fid_signal_range,
        signal_pathway=fid_pathway,
        fl=fl,
        title=f"{fir_node} rough alignment",
    )

    fl.basename = f"{real_file} {ep_node}"
    Ep = psd.find_file(
        real_file,
        exp_type=real_exptype,
        expno=ep_node,
        lookup=prscr.lookup_table,
        fl=fl,
    )
    if Ep.get_prop("log") is None:
        Ep = prscr.attach_log_data_from_file(Ep, real_file, real_exptype)
    m = re.search(r".*ODNP.*v([0-9]+)$", Ep.get_prop("postproc_type"))
    if m is None:
        raise IOError(
            f"Unexpected postproc_type: {Ep.get_prop('postproc_type')!r}"
        )
    if int(m.groups()[0]) < 6:
        Ep = generate_coordinates_from_log(Ep, fl=fl)
    orig_axis_error = Ep.get_error("indirect")
    Ep["indirect"] = Ep["indirect"]["time"]
    Ep.set_units("indirect", "s")
    Ep.set_error("indirect", orig_axis_error["time"])
    ep_pathway = Ep.get_prop("coherence_pathway")
    frq_center, frq_half = prscr.find_peakrange(
        prscr.select_pathway(Ep.C, ep_pathway),
        direct="t2",
        peak_lower_thresh=0.1,
    )
    ep_signal_range = tuple(
        sorted((frq_center - frq_half, frq_center + frq_half))
    )
    prscr.table_of_integrals(
        Ep.C,
        signal_range=ep_signal_range,
        signal_pathway=ep_pathway,
        fl=fl,
        title="E(p) correlation alignment",
        repeat_dims="indirect",
        alignment_sigma=alignment_mask_sigma,
        max_shift=Ep_alignment_max_shift,
        center_aligned_peak=False,
        show_alignment_diagnostics=show_alignment_diagnostics,
        require_unitary_error=False,
        equal_energy_apodization=equal_energy_apodization,
    )
    prscr.rough_table_of_integrals(
        Ep.C,
        signal_range=ep_signal_range,
        signal_pathway=ep_pathway,
        fl=fl,
        title="E(p) rough alignment",
    )
    # }}}
