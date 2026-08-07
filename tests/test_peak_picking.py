import numpy as np
import pytest
import sympy as sp

from pyspecdata import lmfitdata, nddata
from pyspecProcScripts import peel_peaks


def derivative_lineshape(field, center, linewidth, amplitude, balance):
    """Unit-area Lorentzian/Gaussian derivative mixture."""
    offset = field - center
    half_width = linewidth / 2
    lorentzian = (
        -2 * half_width * offset / (np.pi * (offset**2 + half_width**2) ** 2)
    )
    gaussian = (
        np.sqrt(4 * np.log(2) / np.pi)
        / linewidth
        * np.exp(-4 * np.log(2) * offset**2 / linewidth**2)
        * (-8 * np.log(2) * offset / linewidth**2)
    )
    return amplitude * ((1 - balance) * lorentzian + balance * gaussian)


def configured_fit(normalization=1.0, balance_value=0.35):
    field = np.linspace(-30, 30, 6001)
    truth = [(-11, 1.1, 0.8), (12, 1.8, 0.5)]
    spectrum = sum(
        derivative_lineshape(
            field, center, linewidth, amplitude, balance_value
        )
        for center, linewidth, amplitude in truth
    )
    fit = lmfitdata(nddata(spectrum, "B").setaxis("B", field))
    B = sp.symbols("B", real=True)
    amplitudes = sp.symbols("height0:2", real=True)
    linewidths = sp.symbols("width0:2", real=True)
    centers = sp.symbols("center0:2", real=True)
    balances = sp.symbols("balance0:2", real=True)
    expression = 0
    for amplitude, linewidth, center, balance in zip(
        amplitudes, linewidths, centers, balances
    ):
        offset = B - center
        half_width = linewidth / 2
        lorentzian = (
            -2
            * half_width
            * offset
            / (sp.pi * (offset**2 + half_width**2) ** 2)
        )
        gaussian = (
            sp.sqrt(4 * sp.log(2) / sp.pi)
            / linewidth
            * sp.exp(-4 * sp.log(2) * offset**2 / linewidth**2)
            * (-8 * sp.log(2) * offset / linewidth**2)
        )
        expression += (
            normalization
            * amplitude
            * ((1 - balance) * lorentzian + balance * gaussian)
        )
    fit.functional_form = expression
    fit.set_guess({str(balance): balance_value for balance in balances})

    @fit.define_data_transform
    def identity_data_transform(d_local):
        return d_local

    @fit.define_residual_transform
    def identity_residual_transform(d_local):
        return d_local

    return fit, amplitudes, linewidths, centers, balances, truth


@pytest.mark.parametrize("balance_value", [0.1, 0.8])
def test_peel_peaks_calibrates_each_model_line_independently(balance_value):
    fit, amplitudes, linewidths, centers, balances, truth = configured_fit(
        balance_value=balance_value
    )

    returned = peel_peaks(
        fit,
        [str(symbol) for symbol in amplitudes],
        [str(symbol) for symbol in linewidths],
        [str(symbol) for symbol in centers],
        close_threshold=3,
    )

    assert returned is fit
    for j, (center, linewidth, amplitude) in enumerate(truth):
        np.testing.assert_allclose(
            fit.guess_parameters[f"center{j}"].value, center, atol=0.02
        )
        np.testing.assert_allclose(
            fit.guess_parameters[f"width{j}"].value, linewidth, rtol=0.03
        )
        np.testing.assert_allclose(
            fit.guess_parameters[f"height{j}"].value, amplitude, rtol=0.04
        )
        assert fit.guess_parameters[f"height{j}"].min == 0
        assert fit.guess_parameters[f"height{j}"].max > amplitude
        np.testing.assert_allclose(
            fit.guess_parameters[f"width{j}"].min,
            0.1 * fit.guess_parameters[f"width{j}"].value,
        )
        np.testing.assert_allclose(
            fit.guess_parameters[f"width{j}"].max,
            10 * fit.guess_parameters[f"width{j}"].value,
        )
        assert fit.guess_parameters[f"balance{j}"].value == balance_value


def test_peel_peaks_accounts_for_model_normalization():
    fit, amplitudes, linewidths, centers, _, truth = configured_fit(
        normalization=7.0
    )

    returned = peel_peaks(
        fit,
        [str(symbol) for symbol in amplitudes],
        [str(symbol) for symbol in linewidths],
        [str(symbol) for symbol in centers],
        close_threshold=3,
    )

    assert returned is fit
    for j, (_, _, amplitude) in enumerate(truth):
        np.testing.assert_allclose(
            fit.guess_parameters[f"height{j}"].value,
            amplitude / 7.0,
            rtol=0.04,
        )


def test_peel_peaks_returns_empty_dict_for_zero_signal():
    fit, amplitudes, linewidths, centers, _, _ = configured_fit()
    fit *= 0

    returned = peel_peaks(
        fit,
        [str(symbol) for symbol in amplitudes],
        [str(symbol) for symbol in linewidths],
        [str(symbol) for symbol in centers],
        close_threshold=3,
    )

    assert returned is fit
