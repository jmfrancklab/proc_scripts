"""Initial-guess utilities for derivative ESR spectra."""

import numpy as np


def peel_peaks(
    d,
    amplitude_parameters,
    linewidth_parameters,
    center_parameters,
    close_threshold,
    exclusion_factor=2.0,
    lo=0.95,
    hi=0.05,
    n=80,
    axis="B",
):
    r"""Find derivative peaks and install model-calibrated guesses on ``d``.

    ``d`` must be an ``lmfitdata`` with its functional form, data transform,
    residual transform, and non-peak guesses already configured.  Each
    position in ``amplitude_parameters``, ``linewidth_parameters``, and
    ``center_parameters`` identifies the three fit parameters belonging to
    one line.  In particular, ``linewidth_parameters`` names the parameters
    that control overall width, not shape-balance parameters such as a
    Lorentzian/Gaussian fraction.

    The experimental trace is generated with the data transform and peeled
    one line at a time.  Each isolated model line is then evaluated to convert
    its measured peak-to-peak separation into the named linewidth parameter
    and its measured peak-to-peak height into the named amplitude parameter.
    No lineshape-specific conversion factors are used.

    Parameters
    ----------
    d : pyspecdata.lmfitdata
        Configured fit object.  Its guess parameters are updated in place.
    amplitude_parameters : sequence of str
        Amplitude parameter names, ordered by increasing center field.
    linewidth_parameters : sequence of str
        Overall-width parameter names, ordered by increasing center field.
    center_parameters : sequence of str
        Center-field parameter names, ordered by increasing center field.
    close_threshold : float
        Largest accepted positive/negative chunk-center separation, in the
        units of ``axis``.
    exclusion_factor : float, optional
        Half-width of the region masked after finding a peak, in units of its
        measured peak-to-peak separation.
    lo, hi : float, optional
        Starting and ending fractions of the remaining signal maximum used
        for the threshold scan.
    n : int, optional
        Number of thresholds in the scan.
    axis : str, optional
        Name of the field axis.

    Returns
    -------
    pyspecdata.lmfitdata
        The same fit object, with calibrated guesses and bounds installed.
    """
    max_peaks = len(amplitude_parameters)
    experimental = d.data_transform(d.C)
    remaining = experimental.C
    found = []
    for _ in range(max_peaks):
        signal_scale = max(abs(remaining.max()), abs(remaining.min()))
        if signal_scale == 0:
            break
        best_pair = None
        best_distance = None
        for threshold in np.linspace(lo, hi, n):
            cutoff = threshold * signal_scale
            # contiguous returns chunks widest-first.  Zero-width chunks do
            # not contain a lobe, and chunks beyond this small multiple of
            # max_peaks are much more likely to be noise than real lines.
            positive_chunks = remaining.contiguous(
                lambda x, cutoff=cutoff: x > cutoff
            )
            positive_chunks = positive_chunks[
                np.abs(np.diff(positive_chunks, axis=1)).ravel() > 0
            ][: 2 * max_peaks]
            negative_chunks = remaining.contiguous(
                lambda x, cutoff=cutoff: x < -cutoff
            )
            negative_chunks = negative_chunks[
                np.abs(np.diff(negative_chunks, axis=1)).ravel() > 0
            ][: 2 * max_peaks]
            if len(positive_chunks) == 0 or len(negative_chunks) == 0:
                continue
            distances = abs(
                positive_chunks.mean(axis=1)[:, None]
                - negative_chunks.mean(axis=1)[None, :]
            )
            positive_index, negative_index = np.unravel_index(
                np.argmin(distances), distances.shape
            )
            distance = distances[positive_index, negative_index]
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_pair = (
                    positive_chunks[positive_index],
                    negative_chunks[negative_index],
                )
            if best_distance < close_threshold:
                break
        if best_pair is None or best_distance >= close_threshold:
            break

        positive_slice = remaining[axis : best_pair[0]]
        negative_slice = remaining[axis : best_pair[1]]
        positive = positive_slice.argmax()
        positive["y"] = positive_slice.max()
        negative = negative_slice.argmin()
        negative["y"] = negative_slice.min()
        field_extrema = np.array([positive[axis], negative[axis]])
        height_extrema = np.array([positive["y"], negative["y"]])
        peak_to_peak_width = np.ptp(field_extrema)
        center_field = field_extrema.mean()
        found.append(
            {
                "dB_pp": peak_to_peak_width,
                "raw_height": np.ptp(height_extrema),
                "center": center_field,
            }
        )
        exclusion_bounds = center_field + (
            exclusion_factor * peak_to_peak_width * np.array([-1, 1])
        )
        remaining[axis:exclusion_bounds] = 0

    found.sort(key=lambda peak: peak["center"])
    parameter_groups = list(
        zip(amplitude_parameters, linewidth_parameters, center_parameters)
    )
    # Give every line finite parameters before evaluating an isolated
    # one: zero amplitude does not protect a symbolic model from 0 *
    # NaN.
    # TODO ☐: the previous comment is unclear -- where the heck would a NaN come from? Are you considering edge cases that will never occur??!!
    for amplitude_name, linewidth_name, center_name in parameter_groups:
        d.guess_parameters[amplitude_name].value = 0.0
        d.guess_parameters[linewidth_name].value = 1.0
        d.guess_parameters[center_name].value = 0.0
    calibrated = {}
    active_groups = list(zip(found, parameter_groups))
    for peak, (amplitude_name, linewidth_name, center_name) in active_groups:
        d.guess_parameters[amplitude_name].value = 1.0
        d.guess_parameters[linewidth_name].value = peak["dB_pp"]
        d.guess_parameters[center_name].value = peak["center"]
        model_line = d.set_to_guess().eval()
        # TODO ☐: the following is stupid -- just directly subtract the
        #         max from the min.  Also use .max() and .min() rather
        #         than argmax
        model_field_extrema = np.array(
            [model_line.argmax()[axis], model_line.argmin()[axis]]
        )
        model_width = np.ptp(model_field_extrema)
        linewidth = (
            d.guess_parameters[linewidth_name].value
            * peak["dB_pp"]
            / model_width
        )
        d.guess_parameters[linewidth_name].value = linewidth

        model_line = d.set_to_guess().eval()
        # TODO ☐: same comment as above
        model_height = np.ptp([model_line.max(), model_line.min()])
        amplitude = peak["raw_height"] / model_height
        calibrated.update(
            {
                amplitude_name: amplitude,
                linewidth_name: linewidth,
                center_name: peak["center"],
            }
        )
        d.guess_parameters[amplitude_name].value = 0.0

    for name, value in calibrated.items():
        d.guess_parameters[name].value = value
    for _, (amplitude_name, linewidth_name, _) in active_groups:
        amplitude = d.guess_parameters[amplitude_name]
        amplitude.min = 0
        amplitude.max = 10 * max(amplitude.value, 1e-12)
        linewidth = d.guess_parameters[linewidth_name]
        linewidth.min = 0.1 * linewidth.value
        linewidth.max = 10 * linewidth.value
    d.set_to_guess()
    return d
