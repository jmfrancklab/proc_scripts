from pylab import *
from pyspecdata import *
from pyspecProcScripts import Heaviside_time_domain

s = nddata(r_[-5:5:11j], "t")
s.set_ft_prop("t", True)
with figlist_var() as fl:
    # show that we can change by sub-integral amounts
    for thisrange in [(1, 4), (0.75, 3.75)]:
        result = Heaviside_time_domain(s, thisrange, direct="t")
        result.ft("t")
        print(result.C.sum("t").item())
        plot(result, "o-")
