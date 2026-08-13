import numpy as np
import pytest

from pyspecdata import nddata
from pyspecProcScripts import find_peakrange


@pytest.mark.parametrize("data_units", [None, "V"])
def test_find_peakrange_preserves_input_with_optional_units(data_units):
    t2 = np.linspace(-5e-3, 5e-3, 512, endpoint=False)
    data = (
        nddata(np.exp(-(t2 / 1e-3) ** 2), [-1], ["t2"])
        .setaxis("t2", t2)
        .set_units("t2", "s")
    )
    if data_units is not None:
        data.set_units(data_units)
    data.ft("t2", shift=True)
    original = data.C

    center, half_width = find_peakrange(data)

    assert np.isfinite(center)
    assert np.isfinite(half_width)
    assert half_width > 0
    assert data.dimlabels == original.dimlabels
    assert data.get_units() == original.get_units()
    assert data.get_units("t2") == original.get_units("t2")
    np.testing.assert_allclose(data.getaxis("t2"), original.getaxis("t2"))
    np.testing.assert_allclose(data.data, original.data)
