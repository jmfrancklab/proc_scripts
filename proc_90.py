from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from proc_scripts import zeroth_order_ph, calc_baseline, ph1_real_Abs
from proc_scripts.load_data import postproc_dict
from scipy.optimize import minimize,curve_fit,least_squares
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
fl = figlist_var()
for searchstr, exp_type, which_exp, postproc in [
        ('w8_200731','test_equip',1,'zg2h'),
        ]:
    s = find_file(searchstr, exp_type=exp_type,
                expno=which_exp, postproc=postproc, lookup=postproc_dict)     
    fl.show();quit()
    if manual_phcyc:
        fl.basename = exp_name
        fl.next('image of phase cycle domain')
        fl.image(d)
        d.setaxis('ph',r_[0,1,2,3]/4.)
        d.ft('ph')# based on 90 pulse experiment, we seem to want ft rather than ift for deuterium (phcyc goes around the circle the wrong way)
        fl.next('image of coherence domain')
        fl.image(d)
        d = d['ph',-1].C
    fl.basename = None

    baseline = (d['t2':(None,-100)].C.mean('t2', return_error=False)+d['t2':(300,None)].C.mean('t2', return_error=False)).item()*0.5
    d -= baseline

    # Use this to get the 1st order phase correction
    SW = diff(d.getaxis('t2')[r_[0,-1]]).item()
    # {{{ generate a series of test spectra with first order correction applied
    thisph1 = nddata(r_[-6:6:2048j]/SW,'phi1').set_units('phi1','s')
    phase_test_r = d * exp(-1j*2*pi*thisph1*d.fromaxis('t2'))
    # }}}
    # {{{ correct the test spectra by the zeroth order correction
    phase_test_rph0 = phase_test_r.C.sum('t2')
    phase_test_rph0 /= abs(phase_test_rph0)
    phase_test_r /= phase_test_rph0
    # }}}
    # {{{ calculate and plot the cost as a function of the first order correction
    cost = abs(phase_test_r.real).sum('t2')
    fl.next('phasing cost function for first order correction')
    fl.plot(cost,'.',human_units=False)
    # }}}
    # {{{ determine and apply the optimal phase correction
    ph1_opt = cost.argmin('phi1').item()
    d *= exp(-1j*2*pi*ph1_opt*d.fromaxis('t2')) # first order corr
    # {{{ determining 0th order correction
    d_ph0 = d.C.sum('t2') # summing along t2
    d_ph0 /= abs(d_ph0) # dividing by abs val of the sum
    d /= d_ph0
    # }}}
    # }}}
    # {{{ show the corrected, unshifted spectrum
    fl.next('corrected spectrum')
    fl.plot(d,label="%s ph1=%f ms"%(label_id,ph1_opt/1e-3),human_units=False)
    # }}}
    # {{{ normalize, shift the axis, and select the center 300 Hz of the spectrum
    fl.next('normalized and centered')
    max_val_at = d.C.argmax('t2').item() #adding maxes to be put on git screen
    d /= d['t2':max_val_at]
    d.setaxis('t2',lambda x: x-max_val_at)
    print("SFO1",d.get_prop('acq')['SFO2'])
    fl.plot(d['t2':(-150,150)], alpha=0.5, label=r"%s $\nu_0$=%f ppm"%(label_id,max_val_at/d.get_prop('acq')['SFO1']),human_units=False)
    # }}}
fl.show() # showing the plot
quit()
