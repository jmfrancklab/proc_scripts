import numpy as np

from pyspecdata import nddata
from pyspecProcScripts import peel_peaks


def derivative_lorentzian(field, center, fwhm, amplitude):
    """Derivative of a unit-area Lorentzian absorption line."""
    offset = field - center
    half_width = fwhm / 2
    return (
        -2
        * amplitude
        * half_width
        * offset
        / (np.pi * (offset**2 + half_width**2) ** 2)
    )


def test_peel_peaks_finds_sorted_centers_and_pools_linewidth():
    field = np.linspace(-30, 30, 6001)
    spectrum = sum(
        derivative_lorentzian(field, center, 1.2, amplitude)
        for center, amplitude in [(12, 0.6), (-11, 1.0), (1, 0.8)]
    )
    data = nddata(spectrum, "B").setaxis("B", field)
    original = data.C

    peaks = peel_peaks(data, 3, close_threshold=2)

    assert len(peaks) == 3
    np.testing.assert_allclose(
        [peak["Bcenter"] for peak in peaks], [-11, 1, 12], atol=0.02
    )
    assert len({peak["FWHM"] for peak in peaks}) == 1
    assert all(peak["A"] > 0 for peak in peaks)
    np.testing.assert_allclose(data.data, original.data)
    np.testing.assert_allclose(data.getaxis("B"), original.getaxis("B"))


def test_peel_peaks_returns_empty_list_for_zero_signal():
    data = nddata(np.zeros(101), "B").setaxis("B", np.linspace(-5, 5, 101))

    assert peel_peaks(data, 4, close_threshold=1) == []
