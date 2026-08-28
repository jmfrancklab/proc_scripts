import matplotlib.pyplot as plt
import numpy as np
import pytest

import pyspecdata as psd
import pyspecProcScripts.phasing as phasing

from pyspecProcScripts import (
    det_inh_bounds,
    fid_side_from_echo,
    fid_from_echo,
    fit_envelope,
    find_exponential_echo_center,
    zeroth_order_ph,
)


def _synthetic_spectrum(
    peaks=((100.0, 70.0, 1.0),),
    noise=0.0,
    echo_center=None,
    phase=0.0,
    points=2048,
    dwell_time=25e-6,
    seed=1,
):
    """Return a frequency-domain synthetic FID or echo."""
    direct = "t2"
    time_axis = np.arange(points) * dwell_time
    if echo_center is None:
        decay_time = time_axis
    else:
        decay_time = abs(time_axis - echo_center)
    signal = sum(
        amplitude
        * np.exp(
            1j * (2 * np.pi * frequency * time_axis + phase)
            - np.pi * linewidth * decay_time
        )
        for frequency, linewidth, amplitude in peaks
    )
    if noise:
        rng = np.random.default_rng(seed)
        signal += noise * (
            rng.standard_normal(time_axis.size)
            + 1j * rng.standard_normal(time_axis.size)
        )
    return (
        psd.nddata(signal, direct)
        .setaxis(direct, time_axis)
        .set_units("V")
        .set_units(direct, "s")
        .ft(direct, shift=True)
    )


def test_det_inh_bounds_matches_threshold_and_preserves_input():
    direct = "t2"
    center_frequency = 100.0
    linewidth = 20.0
    spectrum = _synthetic_spectrum(peaks=((center_frequency, linewidth, 1.0),))
    original_data = spectrum.data.copy()
    original_axis = spectrum.getaxis(direct).copy()
    original_units = spectrum.get_units()
    original_axis_units = spectrum.get_units(direct)
    original_ft_state = spectrum.get_ft_prop(direct)

    detected_center, detected_half_width = det_inh_bounds(
        spectrum,
        0.1,
        peak_lowest_thresh=0.03,
        smoothing_width=1.0,
        echo_like=False,
    )

    frequency_step = abs(spectrum.get_ft_prop(direct, "df"))
    expected_half_width = linewidth / 2 * np.sqrt(0.03**-2 - 1)
    assert abs(detected_center - center_frequency) <= frequency_step
    assert abs(detected_half_width - expected_half_width) <= frequency_step
    np.testing.assert_allclose(
        spectrum.get_prop("inh_bounds"),
        detected_center + np.array([-1, 1]) * detected_half_width,
    )
    np.testing.assert_allclose(spectrum.data, original_data)
    np.testing.assert_allclose(spectrum.getaxis(direct), original_axis)
    assert spectrum.get_units() == original_units
    assert spectrum.get_units(direct) == original_axis_units
    assert spectrum.get_ft_prop(direct) == original_ft_state


def test_det_inh_bounds_is_stable_with_noise_and_an_isolated_spike():
    clean = _synthetic_spectrum()
    noisy = _synthetic_spectrum(noise=2e-3)
    spiked = clean.C
    spike_index = np.argmin(abs(spiked.getaxis("t2") + 4e3))
    spiked.data[spike_index] += 0.3 * abs(spiked.data).max()

    clean_bounds = np.array(det_inh_bounds(clean, 0.1, echo_like=False))
    frequency_step = abs(clean.get_ft_prop("t2", "df"))
    for spectrum in (noisy, spiked):
        detected_bounds = np.array(
            det_inh_bounds(spectrum, 0.1, echo_like=False)
        )
        np.testing.assert_allclose(
            detected_bounds,
            clean_bounds,
            atol=2 * frequency_step,
        )


def test_det_inh_bounds_merges_nearby_fragments(monkeypatch):
    spectrum = _synthetic_spectrum()
    ranges_by_threshold = iter(
        [
            [np.array([-120.0, -80.0]), np.array([80.0, 120.0])],
            [np.array([-150.0, -60.0]), np.array([60.0, 150.0])],
            [np.array([-180.0, -40.0]), np.array([40.0, 180.0])],
        ]
    )

    def contiguous_ranges(*args, **kwargs):
        return next(ranges_by_threshold)

    monkeypatch.setattr(psd.nddata, "contiguous", contiguous_ranges)

    center, half_width = det_inh_bounds(
        spectrum,
        0.1,
        peak_lowest_thresh=0.03,
        echo_like=False,
    )

    assert center == 0.0
    assert half_width == 180.0


def test_det_inh_bounds_rejects_separate_peaks():
    spectrum = _synthetic_spectrum(
        peaks=((-1500.0, 10.0, 1.0), (1500.0, 10.0, 1.0))
    )

    with pytest.raises(ValueError, match="peak|spectral"):
        det_inh_bounds(spectrum, 0.1, echo_like=False)


def test_det_inh_bounds_compensates_a_stored_digital_filter():
    unfiltered = _synthetic_spectrum()
    filtered = unfiltered.C
    digital_filter = np.exp(1j * 2 * np.pi * filtered.getaxis("t2") * 25e-6)
    filtered.data /= digital_filter
    filtered.set_prop("dig_filter", digital_filter)

    unfiltered_bounds = det_inh_bounds(
        unfiltered,
        0.1,
        echo_like=False,
    )
    filtered_bounds = det_inh_bounds(
        filtered,
        0.1,
        echo_like=False,
    )

    np.testing.assert_allclose(filtered_bounds, unfiltered_bounds)


def test_det_inh_bounds_plotting_does_not_change_bounds():
    without_plotting = _synthetic_spectrum()
    with_plotting = without_plotting.C

    expected = det_inh_bounds(without_plotting, 0.1, echo_like=False)
    figure_list = psd.figlist_var()
    plotted = det_inh_bounds(
        with_plotting,
        0.1,
        echo_like=False,
        fl=figure_list,
    )
    plt.close("all")

    np.testing.assert_allclose(plotted, expected)


def test_det_inh_bounds_echo_uses_consistent_return_signature():
    expected_echo_center = 5.0125e-3
    echo = _synthetic_spectrum(echo_center=expected_echo_center)
    echo_with_plotting = echo.C

    detected = det_inh_bounds(echo, 0.1, echo_like=True, fl=None)
    figure_list = psd.figlist_var()
    plotted = det_inh_bounds(
        echo_with_plotting,
        0.1,
        echo_like=True,
        fl=figure_list,
    )
    plt.close("all")

    assert len(detected) == 2
    assert np.all(np.isfinite(detected))
    assert echo.get_prop("echo_center") is not None
    assert np.isfinite(echo.get_prop("echo_center"))
    assert abs(
        echo.get_prop("echo_center") - expected_echo_center
    ) <= echo.get_ft_prop("t2", "dt")
    np.testing.assert_allclose(plotted, detected)
    np.testing.assert_allclose(
        echo_with_plotting.get_prop("echo_center"),
        echo.get_prop("echo_center"),
    )


@pytest.mark.parametrize(
    "echo_center, linewidth, phase, noise, points",
    [
        (5.0125e-3, 40.0, 0.7, 0.0, 2048),
        (7.3375e-3, 100.0, -1.1, 2e-3, 2048),
        (10.0125e-3, 250.0, 2.0, 5e-3, 4096),
        (4.9875e-3, 500.0, -0.4, 2e-3, 1024),
    ],
)
def test_exponential_echo_center_accuracy(
    echo_center,
    linewidth,
    phase,
    noise,
    points,
):
    echo = _synthetic_spectrum(
        peaks=((100.0, linewidth, 1.0),),
        noise=noise,
        echo_center=echo_center,
        phase=phase,
        points=points,
    )

    detected_center = find_exponential_echo_center(
        echo,
        decay_rate=linewidth,
    )

    assert abs(detected_center - echo_center) <= echo.get_ft_prop("t2", "dt")


def test_exponential_echo_center_preserves_input():
    echo = _synthetic_spectrum(echo_center=5.0125e-3)
    original_data = echo.data.copy()
    original_axis = echo.getaxis("t2").copy()
    original_ft_state = echo.get_ft_prop("t2")

    find_exponential_echo_center(echo)

    np.testing.assert_allclose(echo.data, original_data)
    np.testing.assert_allclose(echo.getaxis("t2"), original_axis)
    assert echo.get_ft_prop("t2") == original_ft_state


@pytest.mark.parametrize(
    "echo, decay_rate, error",
    [
        (
            _synthetic_spectrum(
                peaks=((100.0, 250.0, 1.0),),
                echo_center=50e-3,
            ),
            250.0,
            "boundary",
        ),
        (
            _synthetic_spectrum(echo_center=0.1e-3, points=8),
            40.0,
            "too short",
        ),
    ],
)
def test_exponential_echo_center_rejects_invalid_searches(
    echo,
    decay_rate,
    error,
):
    with pytest.raises(ValueError, match=error):
        find_exponential_echo_center(
            echo,
            decay_rate=decay_rate,
        )


def test_exponential_echo_center_rejects_insufficient_overlap():
    with pytest.raises(ValueError, match="no valid.*search interval"):
        find_exponential_echo_center(
            _synthetic_spectrum(echo_center=5e-3),
            minimum_overlap=1 - 1e-12,
        )


def test_exponential_echo_center_rejects_nonfinite_input():
    echo = _synthetic_spectrum(echo_center=5e-3)
    echo.data[0] = np.nan

    with pytest.raises(ValueError, match="non-finite"):
        find_exponential_echo_center(echo)


def test_fid_side_from_echo_centers_without_altering_the_source():
    direct = "t2"
    echo_center = 5.0125e-3
    echo = _synthetic_spectrum(
        echo_center=echo_center,
        phase=0.7,
    )
    original_data = echo.data.copy()
    original_axis = echo.getaxis(direct).copy()
    original_units = echo.get_units()
    original_axis_units = echo.get_units(direct)
    original_ft_state = echo.get_ft_prop(direct)

    fid_side = fid_side_from_echo(echo, echo_center)

    unweighted_fid_side = echo.C.ift(direct)
    unweighted_fid_side[direct] -= (
        unweighted_fid_side.getaxis(direct)[0] + echo_center
    )
    unweighted_fid_side.register_axis({direct: 0})
    unweighted_fid_side = unweighted_fid_side[direct:(0, None)]
    assert fid_side.getaxis(direct)[0] == 0
    assert not fid_side.get_ft_prop(direct)
    assert fid_side.get_prop("echo_center") == echo_center
    np.testing.assert_allclose(
        fid_side[direct, 0].data * 2,
        unweighted_fid_side[direct, 0].data,
    )
    np.testing.assert_allclose(
        fid_side.data[1:],
        unweighted_fid_side.data[1:],
    )
    np.testing.assert_allclose(echo.data, original_data)
    np.testing.assert_allclose(echo.getaxis(direct), original_axis)
    assert echo.get_units() == original_units
    assert echo.get_units(direct) == original_axis_units
    assert echo.get_ft_prop(direct) == original_ft_state

    round_trip = fid_side.C.ft(direct).ift(direct)
    np.testing.assert_allclose(round_trip.data, fid_side.data)
    np.testing.assert_allclose(
        round_trip.getaxis(direct),
        fid_side.getaxis(direct),
    )


def test_fid_side_from_echo_compensates_stored_digital_filter():
    echo_center = 5.0125e-3
    unfiltered = _synthetic_spectrum(echo_center=echo_center)
    filtered = unfiltered.C
    digital_filter = np.exp(1j * 2 * np.pi * filtered.getaxis("t2") * 25e-6)
    filtered.data /= digital_filter
    filtered.set_prop("dig_filter", digital_filter)

    expected = fid_side_from_echo(unfiltered, echo_center)
    compensated = fid_side_from_echo(filtered, echo_center)

    np.testing.assert_allclose(compensated.data, expected.data)
    np.testing.assert_allclose(
        compensated.getaxis("t2"),
        expected.getaxis("t2"),
    )
    assert compensated.get_prop("dig_filter") is None
    np.testing.assert_allclose(filtered.get_prop("dig_filter"), digital_filter)


def test_det_inh_bounds_reuses_stored_echo_center(monkeypatch):
    echo = _synthetic_spectrum(echo_center=5.0125e-3)
    echo.set_prop("echo_center", 5.0125e-3)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("echo-center detection should not be repeated")

    monkeypatch.setattr(
        phasing,
        "find_exponential_echo_center",
        fail_if_called,
    )

    detected = det_inh_bounds(echo, 0.1, echo_like=True)

    assert np.all(np.isfinite(detected))


@pytest.mark.parametrize("linewidth", [40.0, 80.0, 120.0, 250.0, 500.0])
def test_fit_envelope_recovers_homogeneous_linewidth(linewidth):
    echo_center = 5.0125e-3
    echo = _synthetic_spectrum(
        peaks=((100.0, linewidth, 1.0),),
        noise=2e-3,
        echo_center=echo_center,
        points=4096,
    )

    fitted_linewidth = fit_envelope(
        fid_side_from_echo(echo, echo_center),
        mult_two=True,
    )

    assert fitted_linewidth == pytest.approx(linewidth, rel=0.1)


@pytest.mark.parametrize("frequency_spread", [0.0, 500.0, 1500.0])
def test_homogeneous_fit_is_independent_of_frequency_spread(frequency_spread):
    direct = "t2"
    echo_center = 5.0125e-3
    homogeneous_linewidth = 120.0
    time_axis = np.arange(4096) * 25e-6
    component_frequencies = 100.0 + frequency_spread * np.r_[-1, 0, 1]
    signal = np.exp(
        1j * 2 * np.pi * component_frequencies[:, None] * time_axis[None, :]
        - np.pi * homogeneous_linewidth * abs(time_axis[None, :] - echo_center)
    )
    echo = (
        psd.nddata(signal, ["inh_component", direct])
        .setaxis("inh_component", component_frequencies)
        .setaxis(direct, time_axis)
        .set_units(direct, "s")
        .ft(direct, shift=True)
    )

    fitted_linewidth = fit_envelope(
        fid_side_from_echo(echo, echo_center),
        mult_two=True,
    )

    assert fitted_linewidth == pytest.approx(homogeneous_linewidth, rel=0.1)


def test_fit_envelope_rejects_an_unusable_decay():
    unusable_decay = (
        psd.nddata(np.zeros(20), "t2")
        .setaxis("t2", np.arange(20) * 25e-6)
        .set_units("t2", "s")
    )

    with pytest.raises(ValueError, match="finite nonzero decay"):
        fit_envelope(unusable_decay)


def test_fid_from_echo_stores_homogeneous_linewidth():
    expected_linewidth = 120.0
    processed = fid_from_echo(
        _synthetic_spectrum(
            peaks=((100.0, expected_linewidth, 1.0),),
            noise=2e-3,
            echo_center=5.0125e-3,
            points=4096,
        ),
        {},
    )

    assert processed.get_prop("homogeneous_linewidth") == pytest.approx(
        expected_linewidth,
        rel=0.1,
    )


def test_ordinary_fid_bypasses_echo_center_detection(monkeypatch):
    def fail_if_called(*args, **kwargs):
        raise AssertionError("echo-center detection should not be called")

    monkeypatch.setattr(
        phasing,
        "find_exponential_echo_center",
        fail_if_called,
    )

    detected = det_inh_bounds(
        _synthetic_spectrum(),
        0.1,
        echo_like=False,
    )

    assert np.all(np.isfinite(detected))


def test_zeroth_order_phase_weights_signal_amplitude_by_default():
    expected_phase = 0.37
    signal = np.linspace(1.0, 2.0, 40) * np.exp(1j * expected_phase)
    background_sign = np.tile(np.array([-1, 1]), 5000)
    background = (
        0.1 * background_sign * np.exp(1j * (expected_phase + np.pi / 2))
    )
    data = psd.nddata(np.r_[signal, background], "sample")

    weighted_error = abs(
        np.angle(zeroth_order_ph(data) / np.exp(1j * expected_phase))
    )
    unweighted_error = abs(
        np.angle(
            zeroth_order_ph(data, weighted=False) / np.exp(1j * expected_phase)
        )
    )

    assert weighted_error < 1e-6
    assert weighted_error < unweighted_error
