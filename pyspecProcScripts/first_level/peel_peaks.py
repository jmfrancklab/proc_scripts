"""Initial-guess utilities for derivative ESR spectra."""

import numpy as np


def peel_peaks(
    d,
    max_peaks,
    close_threshold,
    exclusion_factor=2.0,
    lo=0.95,
    hi=0.05,
    n=80,
    axis="B",
):
    r"""Pick derivative-lineshape peaks from ``d`` one at a time.

    At each threshold (scanned from ``lo`` to ``hi``), find the closest
    positive/negative contiguous-chunk pair.  A derivative lineshape has
    positive and negative lobes close together, whereas unrelated noise
    excursions are typically farther apart.  Once a pair is closer than
    ``close_threshold``, accept it and mask ``exclusion_factor`` times its
    peak-to-peak linewidth around its center before looking for another
    peak.

    All returned lines share a pooled linewidth: the mean of the measured
    peak-to-peak separations.  This pooled value is also used to convert
    peak-to-peak heights to amplitudes.  Lines in a hyperfine multiplet
    physically share a linewidth, and the conversion depends on its square,
    so pooling prevents noise in each individual width from being amplified
    into the amplitude guesses.

    Parameters
    ----------
    d : pyspecdata.nddata
        Real, one-dimensional derivative spectrum.  It is not modified.
    max_peaks : int
        Maximum number of peaks to return.  Fewer are returned when no
        remaining chunk pair is closer than ``close_threshold``.
    close_threshold : float
        Largest accepted positive/negative center-to-center separation, in
        the units of ``axis``.
    exclusion_factor : float, optional
        Half-width of the masked region around each found peak, in units of
        that peak's raw peak-to-peak separation.
    lo, hi : float, optional
        Starting and ending fractions of the remaining signal maximum used
        for the threshold scan.
    n : int, optional
        Number of thresholds in the scan.
    axis : str, optional
        Name of the field axis.

    Returns
    -------
    list of dict
        Peak guesses sorted by increasing center field.  Each dictionary has
        the keys ``"A"``, ``"FWHM"``, and ``"Bcenter"``.
    """
    if d.dimlabels != [axis]:
        raise ValueError(
            "peel_peaks requires one-dimensional data whose only axis is "
            f"{axis!r}; got {d.dimlabels!r}"
        )
    if max_peaks < 0 or int(max_peaks) != max_peaks:
        raise ValueError("max_peaks must be a non-negative integer")
    if close_threshold <= 0:
        raise ValueError("close_threshold must be positive")
    if exclusion_factor <= 0:
        raise ValueError("exclusion_factor must be positive")
    if not 0 <= hi <= lo <= 1:
        raise ValueError("lo and hi must satisfy 0 <= hi <= lo <= 1")
    if n < 1 or int(n) != n:
        raise ValueError("n must be a positive integer")

    # A working copy is necessary here: masking each accepted peak is the
    # defining operation of the peeling algorithm, but callers should retain
    # their unmodified spectrum.
    remaining = d.C
    found = []
    for _ in range(int(max_peaks)):
        signal_scale = max(abs(remaining.max()), abs(remaining.min()))
        if signal_scale == 0:
            break
        best_pair = None
        best_distance = None
        for threshold in np.linspace(lo, hi, int(n)):
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
                "Bcenter": center_field,
            }
        )
        exclusion_bounds = center_field + (
            exclusion_factor * peak_to_peak_width * np.array([-1, 1])
        )
        remaining[axis:exclusion_bounds] = 0

    if not found:
        return []

    # Voigt-derivative peak-to-peak conversion, averaged for the 50/50
    # Lorentzian/Gaussian mixture used as the initial lineshape guess:
    #
    # Lorentzian: dB_pp/FWHM = 1/sqrt(3)
    # Gaussian:   dB_pp/FWHM = 1/sqrt(2*ln(2))
    db_pp_per_fwhm = 0.5 * (1 / np.sqrt(3) + 1 / np.sqrt(2 * np.log(2)))
    # For unit integrated amplitude, the corresponding peak-to-peak heights
    # are 3*sqrt(3)/(pi*FWHM**2) and
    # 8*sqrt(2)*ln(2)*exp(-1/2)/(sqrt(pi)*FWHM**2), respectively.
    height_per_amplitude_over_fwhm_squared = 0.5 * (
        3 * np.sqrt(3) / np.pi
        + 8 * np.sqrt(2) * np.log(2) * np.exp(-0.5) / np.sqrt(np.pi)
    )
    common_fwhm = np.mean([peak["dB_pp"] for peak in found]) / db_pp_per_fwhm
    return sorted(
        [
            {
                "A": (
                    peak["raw_height"]
                    * common_fwhm**2
                    / height_per_amplitude_over_fwhm_squared
                ),
                "FWHM": common_fwhm,
                "Bcenter": peak["Bcenter"],
            }
            for peak in found
        ],
        key=lambda peak: peak["Bcenter"],
    )
