import numpy as np

from pyspecdata import nddata
from pyspecProcScripts import cumulant_rms

SAMPLING_COUNTS = (5, 10, 20, 50, 100, 500)
LINEWIDTH = 1.0
TRAVEL = 3.1 * LINEWIDTH
FIELD = np.linspace(-12.0, 12.0, 4001)


def gaussian_derivative(field_axis, center):
    """Gaussian derivative whose width parameter is its standard deviation."""
    offset = (field_axis - center) / LINEWIDTH
    return -offset * np.exp(-0.5 * offset**2) / LINEWIDTH


def spectral_transition(count):
    """Return the shared three-line path sampled at ``count`` endpoints."""
    motion = np.linspace(0.0, 1.0, count)
    spectra = np.array(
        [
            gaussian_derivative(FIELD, -6.0 * LINEWIDTH + TRAVEL * value)
            + gaussian_derivative(FIELD, 0.0)
            + gaussian_derivative(FIELD, 6.0 * LINEWIDTH - TRAVEL * value)
            for value in motion
        ]
    )
    data = nddata(spectra, ["motion", "$B_0$"])
    data.setaxis("motion", motion).setaxis("$B_0$", FIELD)
    return spectra, data


def test_cumulant_rms_converges_with_transition_sampling():
    transitions = {
        count: spectral_transition(count) for count in SAMPLING_COUNTS
    }

    reference_start, reference_end = transitions[SAMPLING_COUNTS[0]][0][
        [0, -1]
    ]
    for spectra, _ in transitions.values():
        np.testing.assert_array_equal(spectra[0], reference_start)
        np.testing.assert_array_equal(spectra[-1], reference_end)

    terminal = {
        count: cumulant_rms(data, "motion").data[-1]
        for count, (_, data) in transitions.items()
    }
    terminal_values = np.array([terminal[count] for count in SAMPLING_COUNTS])
    assert np.all(np.isfinite(terminal_values))

    dense = terminal[500]
    relative_errors = np.abs(terminal_values - dense) / dense
    assert np.all(np.diff(relative_errors) <= 0.0)

    # The five-point chord sum is an accepted coarse-sampling limitation.
    five_spectrum_deficit = (dense - terminal[5]) / dense
    assert 0.05 <= five_spectrum_deficit <= 0.07
    assert abs(terminal[20] - dense) / dense < 0.01
    assert abs(terminal[100] - dense) / dense < 0.001


FIELD_HALF_RANGE = 12.0  # matches FIELD = linspace(-12, 12, ...) above
# points per linewidth: half / base / double / quadruple of a realistic
# minimum acquisition density (nobody acquires below ~5 pts/linewidth)
FIELD_POINTS_PER_LINEWIDTH = (5, 10, 20, 40)
FIXED_MOTION_COUNT = 100  # dense enough that indirect-dim discretization is not a confound


def spectral_transition_field_sampling(field_count, motion_count=FIXED_MOTION_COUNT):
    """Same three-line transition as ``spectral_transition``, but with a
    fixed motion sampling and a variable number of field-axis points."""
    motion = np.linspace(0.0, 1.0, motion_count)
    field = np.linspace(-FIELD_HALF_RANGE, FIELD_HALF_RANGE, field_count)
    spectra = np.array(
        [
            gaussian_derivative(field, -6.0 * LINEWIDTH + TRAVEL * value)
            + gaussian_derivative(field, 0.0)
            + gaussian_derivative(field, 6.0 * LINEWIDTH - TRAVEL * value)
            for value in motion
        ]
    )
    data = nddata(spectra, ["motion", "$B_0$"])
    data.setaxis("motion", motion).setaxis("$B_0$", field)
    return data


def test_cumulant_rms_invariant_to_field_sampling():
    field_counts = [
        int(round(ppl * 2 * FIELD_HALF_RANGE))
        for ppl in FIELD_POINTS_PER_LINEWIDTH
    ]
    terminal = {
        count: cumulant_rms(
            spectral_transition_field_sampling(count), "motion"
        ).data[-1]
        for count in field_counts
    }
    values = np.array([terminal[count] for count in field_counts])
    assert np.all(np.isfinite(values))
    # Halving/doubling/quadrupling the field-axis point count -- all
    # comfortably above the ~5-points-per-linewidth floor of a real
    # acquisition -- changes the terminal cumulant only at the level of
    # floating-point noise (see cumulant_rms_func.py for why this differs
    # from the indirect-dimension convergence test above).
    np.testing.assert_allclose(values, values[0], rtol=1e-8, atol=1e-10)
