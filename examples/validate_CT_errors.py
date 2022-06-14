""" Validate Inactive CT Error
==============================

Estimates the error of the integral of an actual data set of a standard echo
experiment. Three methods of acquiring the error associated with the data are 
compared:

    -   Taking an area along the active coherence transfer (CT) pathway outside of the bandwidth of the signal
        signal and propagating that error to estimate the error associated with the integral.
        (The traditional method of acquiring the error associated with a data set.)
    -   Taking the integral in the inactive CT pathways and propagating to get the error 
        associated with the integral in the active CT.
    -   Taking the standard deviation of many integrals determined by
        integrating over the signal bandwidth of the active CT pathway.
        (Best method when many scans are available)
    
Demonstrates that by propagating the error of the integral in the inactive CTs we still
get a reasonable error within the limits of traditional methods.
"""
from pyspecdata import *
from pylab import *
from matplotlib import *
from pyspecProcScripts import *
from pyspecProcScripts.correlation_alignment import correl_align
import numpy as np
import sympy as s
from collections import OrderedDict
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery
from numpy.random import seed
seed(120)
# sphinx_gallery_thumbnail_number = 4

fl = figlist_var()
signal_pathway = {"ph1": 1, "ph2": 0}
excluded_pathways = [(0, 0), (0, 3)]
colors = ["r", "darkorange", "gold", "g", "c", "b", "m", "lightcoral"]
# {{{ generate the fake data
t2, nScans, ph1, ph2 = s.symbols('t2 nScans ph1 ph2')
data = fake_data(
    (23*(1+1e-8*nScans)*(s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi))),
    OrderedDict([
        ("nScans" , nddata(ones(10000), "nScans")),
        ("ph2" , nddata(r_[0.0,2.0] / 4, "ph2")),
        ("ph1" , nddata(r_[0., 1., 2., 3.] / 4, "ph1")),
        ("t2" , nddata(r_[0:0.085:256j], "t2"))]),
        signal_pathway, scale = 15.0)
# {{{ just have the data phased
data.labels({'ph2':r_[0.0,2.0]/4,'ph1':r_[0.0,1.0,2.0,3.0]/4})
data.reorder(["ph1", "ph2"])
data.ft('t2')
data /= sqrt(ndshape(data)['t2'])*data.get_ft_prop('t2','dt')
data.ift('t2')
# }}}
data = data["t2" : (0,None)]
data["t2":0] *= 0.5
data.ft("t2")
# }}}
# {{{Normalization
int_frq_slice = integrate_limits(select_pathway(data.C, signal_pathway),cutoff=0.1)
s_integral = select_pathway(data['t2':int_frq_slice].C, signal_pathway).integrate('t2')
avg_d = s_integral.C.mean().real.item()
s_integral /= avg_d
data /= avg_d
data.ift("t2")
# }}}
data.ft('t2')
fl.next('Fake Data')
fl.image(data)
error_pathway = (
    set(
        ((j, k) for j in range(ndshape(data)["ph1"]) for k in range(ndshape(data)["ph2"]))
    )
    - set(excluded_pathways)
    - set([(signal_pathway["ph1"], signal_pathway["ph2"])])
)
error_pathway = [{"ph1": j, "ph2": k} for j, k in error_pathway]
# {{{Making lists for all individual inactive pathways to get error
# associated with each one
var_lst = []
std_lst = []
avg_var_lst = []
avg_std_lst = []
#for thispathway in error_pathway:
#    s_thisint, frq_slice_check = integral_w_errors(
#        data.C,
#        signal_pathway,
#        [thispathway],
#        cutoff = 0.1,
#        indirect="nScans",
#        return_frq_slice=True,
#    )
#    assert all(frq_slice_check == int_frq_slice)
#    std = s_thisint.get_error()
#    var = std**2
#    var_lst.append(var)
#    std_lst.append(std)
#    avg_var = var.mean().item()
#    avg_std = sqrt(avg_var)
#    avg_var_lst.append(avg_var)
#    avg_std_lst.append(avg_std)
# }}}
# {{{ Calculating propagated error averaged over all inactive CTs (as the
#     function is meant to be called)
fl.next('test limits')
inact_frq_slice = integrate_limits(
    select_pathway(data.C, signal_pathway),  
    cutoff=0.1, fl=fl)
s = data['t2':inact_frq_slice].C
f = s.getaxis('t2')
df = f[1] - f[0]
errors = []
for j in range(len(error_pathway)):
    s_forerror = select_pathway(s, error_pathway[j])
    if j == 0:
        N2 = ndshape(s_forerror)['t2']
    s_forerror -= s_forerror.C.mean_all_but(['nScans'])
    s_forerror.run(lambda x: abs(x) ** 2 / 2).mean_all_but(['nScans'])
    s_forerror *= df ** 2  # Δf
    s_forerror *= N2
s = select_pathway(s, signal_pathway)
inact_var = s_forerror.data
inact_std = sqrt(s_forerror.data)
avg_inact_var = inact_var.mean().item()
avg_inact_std = inact_std.mean().item()
fl.plot(select_pathway(data.C.mean('nScans'), signal_pathway))
axvline(inact_frq_slice[-1])
axvline(inact_frq_slice[0])
axvline(int_frq_slice[-1])
axvline(int_frq_slice[0])
# }}}
# {{{ Calculating propagated error along active CT on noise slice
act_frq_slice = integrate_limits(
        select_pathway(data.C, signal_pathway),
        cutoff = 0.1, fl=fl)
slice_start = act_frq_slice[-1] + 500 
off_res_frq_slice = (slice_start,None)
s = data['t2' : off_res_frq_slice].C
s = s['t2',0:N2]
f = s.getaxis('t2')
df = f[1] - f[0]
s_forerror = select_pathway(s, signal_pathway)
s_forerror -= s_forerror.C.mean_all_but(['nScans'])
s_forerror.run(lambda x: np.real(x) ** 2).mean_all_but(['nScans'])
s_forerror *= df ** 2
s_forerror *= N2
s = select_pathway(s,signal_pathway)
active_var = s_forerror.data
active_std = sqrt(s_forerror.data)
avg_active_var = active_var.mean().item()
avg_active_std = active_std.mean().item()
axvline(s.getaxis('t2')[0])
axvline(s.getaxis('t2')[-1])
# }}}
# {{{ Calculating the std dev -- error associated with the integrals
expected = s_integral.run(real).mean('nScans',std=True)
expected_std = expected.get_error()
expected_var = expected_std**2
# }}}
# {{{ Plotting Errors
fl.next("Comparison of variances", legend=True)
#for i in range(len(var_lst)):
    #fl.plot(
    #    var_lst[i],
    #    "o",
    #    color=colors[i],
    #    label="on excluded path of %s" % error_pathway[i],
    #)
fl.plot(active_var, "x", 
        label="propagated error from integration of\nactive CT in noise slice")
fl.plot(inact_var, "o", color="brown",
    label="Error from all inactive CTs")
#for i in range(len(std_lst)):
#    axhline(
#        y=avg_var_lst[i],
#        linestyle=":",
#        color=colors[i],
#        label="averaged %s" % error_pathway[i],
#    )
axhline(
    y=avg_active_var,
    linestyle="--",
    label="averaged propagated error\nfrom active CT in noise slice",
)
axhline(
    y=avg_inact_var,
    linestyle="--",
    color="brown",
    label="averaged propagated error\nfrom all inactive CTs",
)
axhline(
    y=expected_var,
    c="k",
    linestyle="-",
    label="traditional std using numpy",
 )
xlabel('nScans')
ylabel('Associated variances')
fl.next('comparison of standard deviations',legend=True)
#for i in range(len(std_lst)):
#    fl.plot(
#        std_lst[i],
#        "o",
#        color=colors[i],
#        label="on excluded path of %s" % error_pathway[i],
#    )
fl.plot(active_std, "x", 
        label="propagated error from integration of\nactive CT in noise slice")
fl.plot(inact_std, "o", color="brown",
    label="Error from all inactive CTs")
#for i in range(len(std_lst)):
#    axhline(
#        y=avg_std_lst[i],
#        linestyle=":",
#        color=colors[i],
#        label="averaged %s" % error_pathway[i],
#    )
axhline(
    y=avg_active_std,
    linestyle="--",
    label="averaged propagated error\nfrom active CT in noise slice",
)
axhline(
    y=avg_inact_std,
    linestyle="--",
    color="orange",
    label="averaged propagated error\nfrom all inactive CTs",
)
axhline(
    y=expected_std,
    c="k",
    linestyle="-",
    label="traditional std using numpy",
)
xlabel('nScans')
ylabel('Associated Std')
#}}}
print("**** *** ***")
print("sample std:",expected_std)
print("**** *** ***")
print("**** *** ***")
print("propagated error along inactive pathways on noise slice (integral_w_errors):",avg_inact_std)
print("**** *** ***")
print("**** *** ***")
print("propagated error along active pathways on noise slice (active_propagation):",avg_active_std)
print("**** *** ***")
plt.axis("tight")
ax = plt.gca()
lims = list(ax.get_ylim())
lims[0] = 0
#ax.set_ylim(lims)
plt.legend()
fl.show()

