import numpy as np
import pyspecdata as psd

from pyspecProcScripts.phasing import find_peakrange, hermitian_function_test


def test_find_peakrange_with_unit_bearing_data():
    direct = "t2"
    time_axis = np.arange(512) * 50e-6
    center_frequency = 100.0
    spectrum = (
        psd.nddata(
            np.exp(
                1j * 2 * np.pi * center_frequency * time_axis
                - np.pi * 70 * time_axis
            ),
            direct,
        )
        .setaxis(direct, time_axis)
        .set_units("V")
        .set_units(direct, "s")
        .ft(direct, shift=True)
    )

    detected_center, detected_half_width = find_peakrange(
        spectrum,
        peak_lower_thresh=0.1,
    )

    assert abs(detected_center - center_frequency) <= spectrum.get_ft_prop(
        direct, "df"
    )
    assert detected_half_width > 0


def test_hermitian_function_uses_axis_for_echo_time():
    direct = "t2"
    dwell_time = 50e-6
    time_axis = np.arange(512) * dwell_time
    echo_center = 5e-3
    echo = (
        psd.nddata(
            np.exp(
                1j * 2 * np.pi * 100 * (time_axis - echo_center)
                - np.pi * 70 * abs(time_axis - echo_center)
            ),
            direct,
        )
        .setaxis(direct, time_axis)
        .set_units("V")
        .set_units(direct, "s")
        .ft(direct, shift=True)
        .ift(direct)
    )

    detected_center = hermitian_function_test(echo, fl=None)

    assert abs(detected_center - echo_center) <= dwell_time
