import numpy as np
from numpy import r_
from pyspecdata import *
from pyspecProcScripts import *
import sympy as s
import matplotlib.pyplot as plt
from collections import OrderedDict
from pylab import *
from cycler import cycler
init_logging(level="debug")
fl = figlist_var()
apo_fn = 'exp'
default_colors = rcParams['axes.prop_cycle'].by_key()['color']
thisls = [':', # from https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
        (0, (3, 1, 1, 1, 1, 1)),
        (0, (5, 1))]
prop_cycle = (cycler(ls=thisls)*3+cycler(color=default_colors[0:9]))()
#{{{Functions
#{{{apodization fn
def exp_apo(s, filter_const):
    """make an exponential function to fit the envelope of the time domain"""
    myfilter = exp(-abs(s.fromaxis('t2'))/filter_const)
    fl.next('apo fn')
    fl.plot(myfilter)
    return myfilter
#}}}
#{{{heaviside hat apodization fn
def heaviside_hat(s, t_slice):
    "make a symmetric sinc fn from the ift of a heaviside hat fn. Here, t_slice defines the area that is equal to 1. The fn retains the axis of the original t2 axis"
    hat = s.fromaxis('t2')
    hat['t2',:] = 0
    hat['t2':t_slice] = 1
    return hat
#}}}
##{{{ generate analytical sinc fn
def integration_sinc(full_s, t_slice, frq_slice, direct = 't2'):
    "make sinc fn that is 1 at t=0 and also 1 in frequency domain over the integration bounds"
    assert not full_s.get_ft_prop(direct), "data needs to be in the time domain!"
    intwidth = frq_slice[1] - frq_slice[0]
    centerfrq = (frq_slice[1] + frq_slice[0])/2
    t = full_s.fromaxis(direct)
    mysinc = t.C.run(lambda x: np.sinc(intwidth*x))
    mysinc *= exp(-1j*2*pi*centerfrq*t)
    mysinc *= intwidth
    mysinc[direct:(None,t_slice[0])] = 0
    mysinc[direct:(t_slice[1],None)] = 0
    return mysinc
#}}}
# {{{ I know how to write a masked mean or std only along 1 dimension, so
#     use numpy apply_along_axis to make it a function that works along 1
#     dimension of multidimensional data
def masked_mean_multi(x, axis=None):
    "Calculates the mean on a 1D axis"
    assert axis is not None
    def masked_mean(x):
        "this only works for 1D data"
        return np.mean(x[isfinite(x)])
    return np.apply_along_axis(masked_mean,axis,x)
def masked_var_multi(x,axis=None):
    "calculates the variance along a 1D axis"
    assert axis is not None
    def masked_var(x):
        "this only works for 1D data"
        if var_has_imag: # take average of variance along real and image
            return np.var(x[isfinite(x)], ddof=1)/2
        else:
            return np.var(x[isfinite(x)], ddof=1)
    return np.apply_along_axis(masked_var,axis,x)
#}}}
# {{{ use a noan mask before applying variance
def f_integration(s, f_range, signal_pathway, ph_cyc = 'ph1', direct = 't2'):
    "integrate the signal in the f domain and set error to the error propagated in the f domain"
    noise =  s.C
    noise[direct:f_range] = nan
    temp = select_pathway(noise,signal_pathway)
    temp.data[:] = nan
    fl.next('freq noise')
    fl.image(noise)
    noise.run(masked_var_multi,'t2')
    noise.run(masked_mean_multi,'ph1')
    justph0noise = s[ph_cyc,0].run(var,direct)/2
    s_int = select_pathway(s,signal_pathway)['t2':(f_range)]
    dν = s_int.get_ft_prop(direct,'df')
    N = ndshape(s_int)[direct]
    s_int.integrate(direct)
    s_int.set_error(sqrt(noise.data * dν**2 * N)) 
    return s_int
# }}}
#{{{ calculate the error propagated in the t domain
def t_error(s, og_data, mult_fn, f_range, signal_pathway, ph_cyc = 'ph1', direct='t2'):
    "variance of original data * dt * integration of f(t)^2"
    dt = np.diff(s.real.C.getaxis(direct)[r_[0,1]]).item()
    integral_mult_fn = abs(mult_fn).C.run(lambda x: x**2).real.integrate(direct)
    assert not og_data.get_ft_prop('t2'), "get in the t domain so I can unitarily transform!"
    og_data.set_ft_prop(direct,'unitary',None)
    og_data.ft('t2',unitary=True)
    t_var = og_data.C
    t_var[direct:f_range] = nan
    temp = select_pathway(t_var, signal_pathway)
    temp.data[:] = nan
    fl.next('masked noise for variance')
    fl.image(t_var)
    t_var.run(masked_var_multi,'t2')
    t_var.run(masked_mean_multi,'ph1')
    t_var /= 4 # theres a difference of a factor of 4 when comparing the causal to the real
    t_err = t_var.data * dt * integral_mult_fn.data
    t_error = sqrt(t_err)
    return t_error
#}}}
#}}}
# {{{ generate fake data
t2, time, repeats, ph1 = s.symbols("t2 time repeats ph1")
signal_pathway = {"ph1": 1}
offset = 100
f_range = tuple(r_[-300, 300] + offset)
excluded_pathways = [(0,2,3)]
n_repeats = 1000 
SNR = 100
data = fake_data(expression = (
            SNR*s.exp(+1j * 2 * s.pi * offset * (t2) - t2 * 36 * s.pi)
            + 1e-10*repeats
            ),
            axis_coords = OrderedDict([
                ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
                ("t2", nddata(r_[0:1.0:1024j], "t2")),
                ("repeats", nddata(r_[0:n_repeats]+1.0, "repeats")),
            ]),
            signal_pathway = {"ph1": 1},
            fake_data_noise_std = 1.0,
            scale=0)
data.reorder(['ph1','repeats','t2'])
data.ft('t2')
#{{{ note that we start at zero, but still need the echolike → causal conversion
data.ift('t2')
data.set_units('t2','s')
data = data['t2':(0,None)]
data *= 2
data['t2':0] *= 0.5 # turning this off causes causal to match and
#                     real/symmetric to mismatch -- I think this is introducing
#                     a correlated error that affects the integral!
#}}}
data.ft('t2')
#{{{ Normalizing data so that integral is 1
s_integral = data['t2':f_range]
s_integral = select_pathway(s_integral,signal_pathway)
avg_d = s_integral.real.integrate('t2').mean('repeats').item()
data /= avg_d
#}}}
# }}}
# {{{ show that f range makes sense, and check for phase noise
fl.next('show dcct')
fl.image(data)
axvline(x=f_range[0], c='w', alpha=0.5)
axvline(x=f_range[1], c='w', alpha=0.5)
# }}}
# {{{ take original f and t range
f_stop_orig = data.getaxis('t2')[-1]
data.ift('t2')
t_stop_orig = data.getaxis('t2')[-1]
# }}}
orig_causal_data = data.C
#{{{ zero fill and make real/symmetric
data.ft('t2',pad = 3*1024)
data.run(np.real)
data.ft_clear_startpoints('t2').set_ft_prop('t2',None)
data.ift('t2',shift = True)
#}}}
#{{{make apodization function and apodized data
apo_data = data.C
fl.next('data with apo fn')
fl.plot(select_pathway(apo_data['repeats',0],signal_pathway))
if apo_fn == 'HH':
    this_apo_fn = heaviside_hat(apo_data.C,(0,0.5))
    apo_data['t2':(0.5,None)] = 0
elif apo_fn == 'exp':    
    this_apo_fn = exp_apo(apo_data,0.025)
    apo_data *= this_apo_fn
fl.next('data with apo fn')
fl.plot(this_apo_fn)
#}}}
apo_data.ft('t2')
data.ft('t2')
#}}}
#I need to use a copy of orig_causal_data because if I don't it just gets thrown back in the loop and yells at me the second time I try to use it
for (thisdata, apo_fn, og_data, var_has_imag,thislabel) in [
        #(data, 1, orig_causal_data.C, False, 'real-symmetric data'),
        (apo_data, this_apo_fn, orig_causal_data.C, False, 'apodized data')
        ]:
   #{{{I am slicing the part of the signal that is nonzero to pass to the frequency-domain integration fn
    withzero_data = thisdata.C
    withzero_data.ift('t2')
    onlynonzero_data = thisdata.C
    onlynonzero_data.ift('t2')
    onlynonzero_data = onlynonzero_data['t2':(-t_stop_orig,t_stop_orig)]
    #}}}
    onlynonzero_data.ft('t2')
    s_int = f_integration(onlynonzero_data.C.real, f_range, signal_pathway)
    #{{{ compare the propagated f error to actual std and avg
    thisavg = np.mean(s_int.data.real)
    thisstd = np.std(s_int.data.real, ddof=1)
    avgofpropagated = sqrt(mean((s_int.get_error())**2))
    avgofpropagated = mean((s_int.get_error()))
    s_int.setaxis('repeats',r_[0:n_repeats])
    fl.next('test error propagation', legend=True)
    thisprop = next(prop_cycle)
    fl.plot(s_int,'o', color=thisprop['color'],label = 'frequency integral %s'%thislabel,
            alpha=0.25)
    axhline(thisavg, color=thisprop['color'], ls='-', label = 'actual avg %s'%thislabel, alpha=0.5)
    axhline(thisavg+thisstd, color=thisprop['color'], ls='-', label = 'actual std %s'%thislabel, alpha=0.5)
    axhline(thisavg-thisstd, color=thisprop['color'], ls='-', alpha=0.5)
    axhline(thisavg+avgofpropagated, color=thisprop['color'],
            ls=thisprop['ls'], label='propagated error in the f domain %s'%thislabel, alpha=0.5)
    axhline(thisavg-avgofpropagated,
            color=thisprop['color'], ls=thisprop['ls'], alpha=0.5)
    #}}}
    onlynonzero_data.ift('t2')
    #{{{ make analytical sinc fn and multiply
    t_lims = onlynonzero_data.C.getaxis('t2')[-1]
    t_range = (-t_lims,t_lims)
    mysinc = integration_sinc(withzero_data, t_slice = t_range, frq_slice = f_range)
    fl.next('mysinc')
    fl.plot(mysinc,label = '%s'%thislabel)
    convolved_data = withzero_data.C * mysinc.C
    #}}}
    #{{{Plot in t domain
    data1d = onlynonzero_data['repeats',0].C
    fl.next('time domain')
    intwidth = f_range[1]-f_range[0]
    fl.plot(select_pathway(data1d,signal_pathway)/abs(select_pathway(data1d.C,signal_pathway).C).max(), 
            label = '%s signal'%thislabel, alpha = 0.5)
    fl.plot(mysinc/abs(mysinc.C).max(), ls = ":", label = '%s sinc fn'%thislabel, alpha = 0.5)
    fl.plot(select_pathway(convolved_data['repeats',0], signal_pathway)/abs(select_pathway(convolved_data['repeats',0],signal_pathway).C).max(),
            label = '%s convolved signal'%thislabel, alpha = 0.5)
    #}}}
    #{{{calculate integral of t domain for convolved signal
    t_int = []
    for j in range(n_repeats):
        thatdata = convolved_data['repeats',j]
        integral_from_t_domain = select_pathway(thatdata.C.real,
                signal_pathway).integrate('t2')
        t_int.append(integral_from_t_domain.real.item())
    t_integral = nddata(np.asarray(t_int), 'repeats')
    t_integral.setaxis('repeats',r_[0:n_repeats])
    #}}}
    #{{{ Plot f domain
    sinc_forplot = mysinc.C
    sinc_forplot.ft('t2')
    data1d.ft('t2')
    convolved_data.ft('t2')
    fl.next('Frequency Domain')
    fl.plot(sinc_forplot*abs(data1d.C).max(), ls = ":", label = '%s sinc fn'%thislabel,
        alpha = 0.5)
    fl.plot(select_pathway(data1d,signal_pathway), label='%s raw data'%thislabel,
            alpha = 0.5)
    fl.plot(select_pathway(convolved_data['repeats',0],signal_pathway)/intwidth,
            label = '%s convolved signal / intwidth'%thislabel, alpha = 0.5)
    plt.axvline(f_range[0],alpha = 0.25)
    plt.axvline(f_range[1],alpha = 0.25)
    #}}}
    #{{{Compare integral of frequency and time domain
    f_int = []
    onlynonzero_data.ft('t2')
    for j in range(n_repeats):
        data = onlynonzero_data['repeats',j]
        integral_from_f_domain = select_pathway(data['t2':f_range].C.real,
                signal_pathway).integrate('t2').item()
        f_int.append(integral_from_f_domain)
    plt.axhline(integral_from_f_domain/intwidth, alpha = 0.25)
    print("congrats! you have proven that the bounded integral in the frequency domain", 
            integral_from_f_domain,  
            "is equal to integral multiplied by the appropriate sinc in the time domain:",
            integral_from_t_domain)
    #}}}
    #{{{ Calculate time domain error
    onlynonzero_data.ift('t2')
    mysinc['t2':0] *= 0.5
    t_errors = t_error(onlynonzero_data, og_data = og_data, mult_fn = mysinc * apo_fn, f_range = f_range, signal_pathway=signal_pathway) #use causal data (original) to calculate the variance
    t_integral.set_error(t_errors)
    t_avg_prop = mean(t_errors)
    #}}}
    #{{{ Plot t integrals and error with f integrals and error
    s_int.setaxis('repeats',r_[0:n_repeats]) #axis is returned where the first point is at 1 instead of 0 so I reset it so everyone is on the same page
    fl.next('test error propagation')
    thisprop = next(prop_cycle)
    fl.plot(t_integral, 'o', color = thisprop['color'], label = 'Time integral %s'%thislabel)
    axhline(thisavg+t_avg_prop, color=thisprop['color'],
            ls=thisprop['ls'], label='propagated error in the t domain', alpha=0.5)
    axhline(thisavg-t_avg_prop, color=thisprop['color'], ls=thisprop['ls'], alpha=0.5)
    ##}}}
fl.show()    



