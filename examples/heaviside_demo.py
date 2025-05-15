""" Heaviside Hat Function Demo
===============================
This example demonstrates the use of the `heaviside_time_domain` function
which generates a Heaviside hat function with a width equal to the range
given,
and then returns the IFT of the Heaviside hat
(for use in time-domain integration).

It specifically deals with frequency slices where the bounds
land somewhere in the middle of the
:math:`\\Delta \\nu`
frequency slices that "belong" to each
discrete point,
and does so by by scaling the points at the
(:math:`\\Delta \\nu_I`)
slice boundary to some fraction of 1
(the fraction to which the integral bounds
enter the frequency slice belonging to the
discrete point).
"""

from pylab import legend
from pyspecdata import nddata, figlist_var, r_, plot
from pyspecProcScripts import heaviside_time_domain

s = nddata(r_[-5:5:11j], "t")
s.set_ft_prop("t", True)
with figlist_var() as fl:
    # show that we can change by sub-integral amounts
    for thisrange in [(0.5, 4.5), (0, 4), (-0.25, 3.75)]:
        result = heaviside_time_domain(s, thisrange, direct="t")
        result.ft("t")
        slice_label = ",".join(str(val) for val in thisrange)
        plot(
            result,
            "o-",
            label="ν slice fed to function: (%s) \n integral: %.3f"
            % (slice_label, result.C.integrate("t").real.item()),
            alpha=0.4,
        )
        legend()
