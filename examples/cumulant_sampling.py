r"""Sampling convergence of the cumulant RMS
===========================================

The cumulant RMS is a dimensionless path-dependent diagnostic.  For a
sequence of spectra it adds the distances between adjacent normalized
spectra, so a finite data set measures the path using straight chords.  As
the same spectral transition is sampled more densely, that discrete sum
converges toward the continuous path length.

Here, two outer Gaussian-derivative lines move inward while a central line
remains fixed.  The ``linewidth`` below is the Gaussian standard deviation
:math:`\sigma`; each moving line travels :math:`3.1\sigma`.  Thus its final
Gaussian envelope at its original center is only
:math:`\exp(-3.1^2/2)\simeq 0.008`, making the initial and final profiles
effectively non-overlapping.  Every sampling count includes both endpoints
of exactly the same motion.

The five-spectrum estimate is deliberately very coarse and underestimates
the dense result by about 5--6 percent.  This is a discretization effect, not
a change in physical units with the number of samples.
"""

import matplotlib.pyplot as plt
import numpy as np
from pyspecdata import nddata

from pyspecProcScripts import cumulant_rms

sampling_counts = (5, 10, 20, 50, 100, 500)
linewidth = 1.0
travel = 3.1 * linewidth
field = np.linspace(-12.0, 12.0, 4001)


def gaussian_derivative(field_axis, center):
    """Return a Gaussian derivative with standard deviation ``linewidth``."""
    offset = (field_axis - center) / linewidth
    return -offset * np.exp(-0.5 * offset**2) / linewidth


def spectral_transition(count):
    """Generate ``count`` spectra along the shared three-line transition."""
    motion = np.linspace(0.0, 1.0, count)
    spectra = np.array(
        [
            gaussian_derivative(field, -6.0 * linewidth + travel * value)
            + gaussian_derivative(field, 0.0)
            + gaussian_derivative(field, 6.0 * linewidth - travel * value)
            for value in motion
        ]
    )
    data = nddata(spectra, ["motion", "$B_0$"])
    data.setaxis("motion", motion).setaxis("$B_0$", field)
    return motion, spectra, data


transitions = {count: spectral_transition(count) for count in sampling_counts}
cumulants = {
    count: cumulant_rms(data, "motion")
    for count, (_, _, data) in transitions.items()
}

# The function returns one value per transition.  For plotting only, put each
# cumulative value at the upper endpoint of its transition and prepend the
# starting point (motion=0, cumulant=0).
plot_cumulants = {
    count: np.concatenate(([0.0], cumulants[count].data))
    for count in sampling_counts
}

figure, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)

coarse_motion, coarse_spectra, _ = transitions[5]
motion_colors = plt.get_cmap("plasma")(
    np.linspace(0.05, 0.9, coarse_motion.size)
)
for value, spectrum, color in zip(
    coarse_motion, coarse_spectra, motion_colors
):
    axes[0].plot(
        field / linewidth,
        spectrum,
        color=color,
        label=rf"$u={value:.2f}$",
    )
axes[0].set(
    xlabel=r"field / Gaussian linewidth $\sigma$",
    ylabel="derivative signal (a.u.)",
    title="Five-spectra transition",
)
# TODO ☐:  do NOT mess with font or marker or linewidth sizes, in general.  Complicates code and make appearance nonstandard!
axes[0].legend(fontsize="small")

count_colors = plt.get_cmap("viridis")(
    np.linspace(0.05, 0.9, len(sampling_counts))
)
for count, color in zip(sampling_counts, count_colors):
    motion = transitions[count][0]
    axes[1].plot(
        motion,
        plot_cumulants[count],
        "o-",
        alpha=0.5,
        # TODO ☐:  do not use a manual color cycle -- stay with matoplotlab default.  Match bar color to line color using the get_plot_color, etc, from pyspec!
        color=color,
        label=f"{count}",
    )
axes[1].set(
    xlabel="motion parameter $u$",
    ylabel="cumulant RMS",
    title="Convergence along the path",
)
axes[1].legend(title="spectra", fontsize="small")

terminal_cumulants = [plot_cumulants[count][-1] for count in sampling_counts]
axes[2].bar(
    [str(count) for count in sampling_counts],
    terminal_cumulants,
    color=count_colors,
)
axes[2].set(
    xlabel="number of spectra",
    ylabel="terminal cumulant RMS",
    title="Convergence of path length",
)

plt.show()
