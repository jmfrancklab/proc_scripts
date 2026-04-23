import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

psd = pytest.importorskip("pyspecdata")

import pyspecProcScripts as prscr


def test_generate_coordinates_from_log_updates_indirect_axis(monkeypatch):
    monkeypatch.setattr(plt, "axvspan", lambda *args, **kwargs: None)

    rng = np.random.default_rng(20260423)
    random_field_axis = rng.uniform(3400.0, 3600.0, size=5)
    field_trace = psd.nddata(random_field_axis.copy(), "field").setaxis(
        "field", random_field_axis
    )
    field_mean = field_trace.mean("field", std=True)
    print("random field axis before mean:", field_trace.getaxis("field"))
    print("field value after mean(std=True):", field_mean.data)
    print("field error after mean(std=True):", field_mean.get_error())

    indirect_dtype = np.dtype(
        [("start_times", np.float64), ("stop_times", np.float64)]
    )
    indirect_axis = np.array(
        [(10.0, 12.0), (12.0, 14.0)],
        dtype=indirect_dtype,
    )
    log_dtype = np.dtype(
        [
            ("time", np.float64),
            ("Rx", np.float64),
            ("power", np.float64),
            ("cmd", "U8"),
        ]
    )
    log = np.array(
        [
            (10.0, 0.0, 0.0, "set"),
            (11.0, 0.0, 0.0, "set"),
            (12.0, 0.0, 0.0, "set"),
            (13.0, 0.0, 0.0, "set"),
        ],
        dtype=log_dtype,
    )

    s = psd.nddata(np.array([1.0, 2.0]), "indirect")
    s.setaxis("indirect", indirect_axis)
    s.set_prop("log", log.copy())

    returned = prscr.generate_coordinates_from_log(
        s,
        filename="unused.h5",
        exp_type="unused",
    )

    assert returned is s
    expected_power = prscr.dBm2power(22.0)
    np.testing.assert_allclose(
        np.asarray(returned.getaxis("indirect")),
        np.array([expected_power, expected_power]),
    )
    np.testing.assert_allclose(
        np.asarray(returned.get_error()),
        np.zeros(2),
    )
