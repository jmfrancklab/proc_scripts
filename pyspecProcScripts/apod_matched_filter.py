from pyspecdata import *


def apod_matched_filter(
    s, axis="t2", convolve_method="gaussian", ret_width=False, fl=None
):
    r"""
    Apodization Matched Filter
    ============================

    This function finds a
    Gaussian or Lorentzian lineshape whose
    linewidth matches that of the data, then
    multiplies the data by this matched-width
    lineshape, ultimately returning an
    apodized version of the input data .

    axis: str
        find matched filter along `axisname`

    convolve_method: str
        specify as one of the following 2 options:
        * Option 1 - 'gaussian'
        finds a matched-filter Gaussian to the
        input data
        * Option 2 - 'lorentzian'
        finds a matched-filter Lorentzian to the
        input data

    ret_width: float
        returns the linewidth found for the matched filter, typically not necessary

    fl: None or figlist_var()
        to show diagnostic plots, set `fl` to the
        figure list; set `fl` to None in order not
        to see any diagnostic plots
    """

    temp = s.C
    sigma = nddata(np.linspace(1e-10, 1e-1, 1000), "sigma").set_units(
        "sigma", "s"
    )
    if convolve_method == "gaussian":
        convolution_set = np.exp(-temp.fromaxis(axis) ** 2 / 2 / sigma**2)
    elif convolve_method == "lorentzian":
        convolution_set = np.exp(-abs(temp.fromaxis(axis)) / sigma)
    signal_E = (abs(temp * convolution_set) ** 2).sum(axis)
    signal_E /= abs(signal_E.data).max()
    if convolve_method == "gaussian":
        filter_width = abs(signal_E - 1 / sqrt(2)).argmin("sigma").item()
    elif convolve_method == "lorentzian":
        filter_width = abs(signal_E - 1 / 2).argmin("sigma").item()
    logger.debug(strm("FILTER WIDTH IS", filter_width))
    if fl is not None:
        fl.next("matched filter diagnostic -- signal Energy")
        fl.plot(signal_E, human_units=False)
        fl.plot(
            signal_E["sigma" : (filter_width, filter_width + 1e-6)],
            "o",
            human_units=False,
        )
    if convolve_method == "gaussian":
        filter_data = np.exp(-temp.fromaxis(axis) ** 2 / 2 / filter_width**2)
    elif convolve_method in "lorentzian":
        filter_data = np.exp(-abs(temp.fromaxis(axis)) / filter_width)
    if fl is not None:
        fl.next("matched filter diagnostic -- time domain")
        fl.plot(abs(temp), alpha=0.6, label="abs val, before mult")
        norm_cnst = (
            abs(temp)[axis : (0 + temp.get_ft_prop(axis, "dt"))]
            / abs(filter_data)[axis : (0 + temp.get_ft_prop(axis, "dt"))]
        )
        fl.plot(
            abs(filter_data) * norm_cnst,
            alpha=0.6,
            label="matched apod. filter",
        )
    temp *= filter_data
    if fl is not None:
        fl.next("matched filter diagnostic -- time domain")
        fl.plot(abs(temp), alpha=0.6, label="abs val, after mult")
    if ret_width:
        return temp, filter_width
    else:
        return temp
