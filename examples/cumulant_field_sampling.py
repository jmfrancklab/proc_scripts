r"""Field-axis sampling invariance of the cumulant RMS
=======================================================

The cumulant RMS diagnostic normalizes each spectrum by its own RMS
amplitude and then takes an RMS difference between adjacent spectra --
both computed as a plain sample mean over the field axis (see
:ref:`sphx_glr_auto_examples_cumulant_sampling.py` for the companion
example, which instead varies the number of spectra along the *indirect*
axis).

A sample-mean quadrature of a smooth lineshape that decays to zero well
inside the field range converges *exponentially* in the number of points,
once the linewidth is resolved.  In practice this means that, unlike the
indirect axis (where the cumulant is a chord approximation to a path length
and only converges slowly, from below, as more spectra are added), the
field axis needs only a modest sampling density -- comfortably below what
any real acquisition would use -- before the cumulant becomes invariant to
the number of field-axis points at essentially floating-point precision.

The left panel shows what varying the field-axis point count actually looks
like: the same single spectrum (at motion parameter :math:`u=0`) sampled at
each field-axis density, plotted as unconnected points so the different
densities are visible on top of one another.  The middle panel instead
varies the *indirect* (motion) coordinate :math:`u` -- the dimension the
cumulant is actually accumulated along -- at a single, well-resolved field
sampling, to show what that variation looks like.  The right panel shows
that the terminal cumulant is the same regardless of which field-axis
density was used.
"""

import matplotlib.pyplot as plt
import numpy as np
from pyspecdata import nddata

from pyspecProcScripts import cumulant_rms

field_half_range = 12.0
points_per_linewidth = (5, 10, 20, 40)
linewidth = 1.0
travel = 3.1 * linewidth
motion_count = 100


def gaussian_derivative(field_axis, center):
    """Return a Gaussian derivative with standard deviation ``linewidth``."""
    offset = (field_axis - center) / linewidth
    return -offset * np.exp(-0.5 * offset**2) / linewidth


def spectrum_at(field, value):
    """Evaluate the shared three-line transition at motion parameter ``value``."""
    return (
        gaussian_derivative(field, -6.0 * linewidth + travel * value)
        + gaussian_derivative(field, 0.0)
        + gaussian_derivative(field, 6.0 * linewidth - travel * value)
    )


def spectral_transition(field_count):
    """Generate the shared three-line transition at a fixed, dense motion
    sampling, with ``field_count`` points along the field axis."""
    motion = np.linspace(0.0, 1.0, motion_count)
    field = np.linspace(-field_half_range, field_half_range, field_count)
    spectra = np.array([spectrum_at(field, value) for value in motion])
    data = nddata(spectra, ["motion", "$B_0$"])
    data.setaxis("motion", motion).setaxis("$B_0$", field)
    return field, data


field_counts = [
    int(round(ppl * 2 * field_half_range)) for ppl in points_per_linewidth
]
transitions = {count: spectral_transition(count) for count in field_counts}
cumulants = {
    count: cumulant_rms(data, "motion")
    for count, (_, data) in transitions.items()
}

figure, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)

sampling_colors = plt.get_cmap("plasma")(
    np.linspace(0.05, 0.9, len(field_counts))
)
for count, ppl, color in zip(
    field_counts, points_per_linewidth, sampling_colors
):
    field = transitions[count][0]
    axes[0].plot(
        field / linewidth,
        spectrum_at(field, 0.0),
        ".",
        alpha=0.5,
        color=color,
        label=f"{ppl} pts/linewidth ({count} field pts)",
    )
axes[0].set(
    xlabel=r"field / Gaussian linewidth $\sigma$",
    ylabel="derivative signal (a.u.)",
    title="Field-axis point density\n(single spectrum, $u=0$)",
)
axes[0].legend()

finest_count = field_counts[-1]
finest_field = transitions[finest_count][0]
demo_motion_values = (0.0, 0.5, 1.0)
demo_linestyles = ("--", ":", "-.")
for value, style in zip(demo_motion_values, demo_linestyles):
    axes[1].plot(
        finest_field / linewidth,
        spectrum_at(finest_field, value),
        style,
        color="C0",
        label=rf"$u={value:.2f}$",
    )
axes[1].set(
    xlabel=r"field / Gaussian linewidth $\sigma$",
    ylabel="derivative signal (a.u.)",
    title="cumulants calculated along\ndimension indicated by $u$",
)
axes[1].legend()

terminal_cumulants = [cumulants[count].data[-1] for count in field_counts]
axes[2].bar(
    [f"{ppl}" for ppl in points_per_linewidth],
    terminal_cumulants,
    color=sampling_colors,
)
axes[2].set(
    xlabel="field points per linewidth",
    ylabel="terminal cumulant RMS",
    title="Invariant to field-axis sampling",
)

plt.show()
