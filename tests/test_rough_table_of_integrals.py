import importlib

import matplotlib.pyplot as plt
import numpy as np
import pytest

import pyspecdata as psd

from pyspecProcScripts import rough_table_of_integrals

table_module = importlib.import_module(
    "pyspecProcScripts.third_level.rough_table_of_integrals"
)


def _synthetic_shifted_fids(frequencies=(-60.0, 40.0)):
    """Return repeated FIDs with known frequency offsets."""
    direct = "t2"
    time_axis = np.arange(1024) * 25e-6
    data = np.exp(
        1j * 2 * np.pi * np.asarray(frequencies)[:, None] * time_axis[None, :]
        - np.pi * 60.0 * time_axis[None, :]
    )
    return (
        psd.nddata(data, ["repeat", direct])
        .setaxis("repeat", np.arange(len(frequencies)))
        .setaxis(direct, time_axis)
        .set_units("repeat", "s")
        .set_units(direct, "s")
        .set_prop("coherence_pathway", {})
        .ft(direct, shift=True)
    )


def _isolate_bounds_flow(monkeypatch):
    """Remove phase/sign changes that are unrelated to bounds migration."""

    def identity_sign(data, *args, **kwargs):
        return data.C.sum("t2").run(lambda x: np.ones_like(x))

    monkeypatch.setattr(table_module, "zeroth_order_ph", lambda data: 1)
    monkeypatch.setattr(table_module, "determine_sign", identity_sign)


@pytest.mark.parametrize("with_plotting", [False, True])
def test_echo_processing_keeps_full_bandwidth_until_fid_from_echo(
    monkeypatch,
    with_plotting,
):
    _isolate_bounds_flow(monkeypatch)
    data = _synthetic_shifted_fids()
    direct = "t2"
    original_axis = data.getaxis(direct).copy()
    frequency_step = abs(data.get_ft_prop(direct, "df"))
    inhomogeneous_bounds = np.array([-250.0, 250.0])
    captured = {}

    def fail_if_called(*args, **kwargs):
        raise AssertionError("supplied bounds should not be redetected")

    def capture_fid_input(fid_data, signal_pathway, **kwargs):
        captured["fid_axis"] = fid_data.getaxis(direct).copy()
        captured["signal_pathway"] = signal_pathway
        captured.update(kwargs)
        return fid_data

    original_integrate = psd.nddata.integrate

    def capture_integration_axis(integral_data, axis, *args, **kwargs):
        captured["integration_axis"] = integral_data.getaxis(axis).copy()
        return original_integrate(integral_data, axis, *args, **kwargs)

    monkeypatch.setattr(table_module, "det_inh_bounds", fail_if_called)
    monkeypatch.setattr(table_module, "fid_from_echo", capture_fid_input)
    monkeypatch.setattr(psd.nddata, "integrate", capture_integration_axis)
    figure_list = psd.figlist_var() if with_plotting else None

    integrals, final_axis = rough_table_of_integrals(
        data,
        signal_range=inhomogeneous_bounds,
        signal_pathway={},
        fl=figure_list,
    )
    plt.close("all")

    np.testing.assert_allclose(captured["fid_axis"], original_axis)
    np.testing.assert_allclose(
        captured["inh_bounds"],
        inhomogeneous_bounds,
    )
    assert captured["signal_pathway"] == {}
    assert captured["max_alignment_shift"] == pytest.approx(
        60.0,
        abs=frequency_step,
    )
    assert captured["integration_axis"][0] >= inhomogeneous_bounds[0]
    assert captured["integration_axis"][-1] <= inhomogeneous_bounds[1]
    assert direct not in integrals.dimlabels
    assert (final_axis is not None) == with_plotting


def test_rough_integrals_detect_bounds_once(monkeypatch):
    _isolate_bounds_flow(monkeypatch)
    data = _synthetic_shifted_fids()
    detected_bounds = np.array([-250.0, 250.0])
    calls = {"detector": 0, "fid": 0}

    def detect_once(detector_data, *args, **kwargs):
        calls["detector"] += 1
        assert kwargs["echo_like"]
        assert kwargs["signal_pathway"] == {}
        detector_data.set_prop("inh_bounds", detected_bounds)
        return 0.0, 250.0

    def reuse_bounds(fid_data, signal_pathway, **kwargs):
        calls["fid"] += 1
        np.testing.assert_allclose(kwargs["inh_bounds"], detected_bounds)
        return fid_data

    monkeypatch.setattr(table_module, "det_inh_bounds", detect_once)
    monkeypatch.setattr(table_module, "fid_from_echo", reuse_bounds)

    rough_table_of_integrals(data, signal_pathway={}, fl=None)

    assert calls == {"detector": 1, "fid": 1}


def test_ordinary_fid_bypasses_echo_processing(monkeypatch):
    _isolate_bounds_flow(monkeypatch)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("ordinary FIDs must bypass fid_from_echo")

    monkeypatch.setattr(table_module, "fid_from_echo", fail_if_called)

    integrals, final_axis = rough_table_of_integrals(
        _synthetic_shifted_fids(),
        signal_range=(-250.0, 250.0),
        signal_pathway={},
        echo_like=False,
        fl=None,
    )

    assert "t2" not in integrals.dimlabels
    assert final_axis is None
