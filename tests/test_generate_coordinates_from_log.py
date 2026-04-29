#!/usr/bin/env python3
# TODO: Delete this script before merging.
# This script tests the ODNP v6-specific log handling path using the
# cleaner_log_save-style serialized log property. The delegated
# proc_spincore_generalproc_v1 processing is stubbed out so this script does
# not test the shared v1 processing path.
import sys
import matplotlib

if __name__ != "__main__" or "--no-show" in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

psd = pytest.importorskip("pyspecdata")

import pyspecProcScripts as prscr
from pyspecProcScripts import load_data as load_data_mod


class SimpleFigList:
    def __init__(self):
        self.fig = None
        self.ax = None

    def next(self, title):
        self.fig, self.ax = plt.subplots(figsize=(8, 4.5))
        self.ax.set_title(title)
        self.ax.set_xlabel("time / s")
        self.ax.set_ylabel("power / W")

    def plot(self, d, fmt="-", human_units=False):
        if self.ax is None:
            self.next("generate_coordinates_from_log test")
        x = np.asarray(d.getaxis(d.dimlabels[0]))
        y = np.asarray(d.data)
        err = d.get_error()
        if err is None:
            self.ax.plot(x, y, fmt)
        else:
            self.ax.errorbar(x, y, yerr=np.asarray(err), fmt=fmt, capsize=4)


def build_fake_dataset():
    indirect_dtype = np.dtype(
        [("start_times", np.float64), ("stop_times", np.float64)]
    )
    indirect_axis = np.array(
        [(10.0, 12.0), (12.0, 14.0)],
        dtype=indirect_dtype,
    )
    expected_log_dtype = np.dtype(
        [
            ("time", np.float64),
            ("Rx", np.float64),
            ("power", np.float64),
            ("field", np.float64),
            ("cmd", np.int64),
        ]
    )
    expected_log = np.array(
        [
            (10.0, 0.0, 0.0, 3400.0, 0),
            (11.0, 0.0, 0.0, 3400.0, 0),
            (12.0, 0.0, 0.0, 3400.0, 0),
            (13.0, 0.0, 0.0, 3400.0, 0),
        ],
        dtype=expected_log_dtype,
    )
    serialized_log_state = {
        "NUMPY_DATA": expected_log,
        "dictkeys": [0],
        "dictvalues": [""],
    }

    s = psd.nddata(np.array([1.0, 2.0]), "indirect")
    s.setaxis("indirect", indirect_axis)
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_prop("log", serialized_log_state)
    return s, expected_log


def run_generate_coordinates_from_log(monkeypatch=None, fl=None):
    if monkeypatch is not None:
        monkeypatch.setattr(plt, "axvspan", lambda *args, **kwargs: None)
        monkeypatch.setattr(
            load_data_mod,
            "proc_spincore_generalproc_v1",
            lambda s, fl=None: s,
        )
    else:
        original = load_data_mod.proc_spincore_generalproc_v1
        load_data_mod.proc_spincore_generalproc_v1 = lambda s, fl=None: s
    s, expected_log = build_fake_dataset()
    try:
        processed = load_data_mod.proc_spincore_ODNP_v6(s)
        assert processed is s
        assert hasattr(processed.get_prop("log"), "total_log")
        np.testing.assert_array_equal(
            np.array(processed.get_prop("log").total_log),
            np.array(
                expected_log,
                dtype=processed.get_prop("log").total_log.dtype,
            ),
        )
        returned = prscr.generate_coordinates_from_log(
            processed,
            filename="unused.h5",
            exp_type="unused",
            fl=fl,
        )
    finally:
        if monkeypatch is None:
            load_data_mod.proc_spincore_generalproc_v1 = original
    return returned


def test_generate_coordinates_from_log_with_proc_spincore_odnp_v6(
    monkeypatch,
):
    returned = run_generate_coordinates_from_log(monkeypatch=monkeypatch)
    assert returned is not None
    expected_power = prscr.dBm2power(22.0)
    np.testing.assert_allclose(
        np.asarray(returned.getaxis("indirect")),
        np.array([expected_power, expected_power]),
    )
    np.testing.assert_allclose(
        np.asarray(returned.get_error()),
        np.zeros(2),
    )


def main():
    fl = SimpleFigList()
    returned = run_generate_coordinates_from_log(fl=fl)
    expected_power = prscr.dBm2power(22.0)
    np.testing.assert_allclose(
        np.asarray(returned.getaxis("indirect")),
        np.array([expected_power, expected_power]),
    )
    np.testing.assert_allclose(
        np.asarray(returned.get_error()),
        np.zeros(2),
    )
    fl.ax.grid(alpha=0.2)
    print(
        "test_generate_coordinates_from_log_with_proc_spincore_odnp_v6 passed"
    )
    if "--no-show" not in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
