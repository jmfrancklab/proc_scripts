""" Heaviside Hat Function Demo
===============================
This example demonstrates the use of the Heaviside_time_domain function
which generates a Heaviside hat function with a width equal to the
range given. Specifically, ranges can be fed both as integral values
or nonintegral values.
"""
from pylab import legend
from pyspecdata import nddata, figlist_var, r_, plot
from pyspecProcScripts import Heaviside_time_domain

s = nddata(r_[-5:5:11j], "t")
s.set_ft_prop("t", True)
with figlist_var() as fl:
    # show that we can change by sub-integral amounts
    for thisrange in [(0, 4), (-0.75, 3.75)]:
        result = Heaviside_time_domain(s, thisrange, direct="t")
        result.ft("t")
        slice_label = ",".join(str(val) for val in thisrange)
        plot(
            result,
            "o-",
            label="ν slice fed to function: (%s) \n sum: %.3f"
            % (slice_label, result.C.sum("t").real.item()),
        )
        legend()
