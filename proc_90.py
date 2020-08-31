from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from proc_scripts import zeroth_order_ph, calc_baseline, ph1_real_Abs
from proc_scripts.load_data import postproc_dict
from scipy.optimize import minimize,curve_fit,least_squares
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
fl = figlist_var()
for searchstr, exp_type, which_exp, postproc, manual_phcyc, fl in [
        ('w8_200731','NMR_Data_AG',1,'zg2h',False,fl),
        ]:
    s = find_file(searchstr, exp_type=exp_type,
                expno=which_exp, postproc=postproc, lookup=postproc_dict,fl=fl)     
    if manual_phcyc:
        fl.basename = exp_name
        fl.next('image of phase cycle domain')
        fl.image(s)
        s.setaxis('ph',r_[0,1,2,3]/4.)
        s.ft('ph')# based on 90 pulse experiment, we seem to want ft rather than ift for deuterium (phcyc goes around the circle the wrong way)
        fl.next('image of coherence domain')
        fl.image(s)
        s = s['ph',-1].C
    #fl.show();quit()
    s.ft('t2')
    baseline = (s['t2':(None,-30)].C.mean('t2')+s['t2':(30,None)].C.mean('t2')).item()*0.5
    s -= baseline
    fl.next('baseline subtracted')
    fl.plot(s)
    #fl.show();quit()
    # Use this to get the 1st order phase correction
    fl.next('first order phase correction')
    dw = diff(s.getaxis('t2')[r_[0,-1]]).item()
    s = ph1_real_Abs(s,dw,fl=fl)
    fl.next('first order applied')
    fl.plot(s)
    # }}}
    # {{{ normalize, shift the axis, and select the center 300 Hz of the spectrum
    fl.next('normalized and centered')
    max_val_at = s.C.argmax('t2').item() #adding maxes to be put on git screen
    s /= s['t2':max_val_at]
    s.setaxis('t2',lambda x: x-max_val_at)
    print("SFO1",s.get_prop('acq')['SFO2'])
    fl.plot(s['t2':(-150,150)], alpha=0.5, label=r"%s $\nu_0$=%f ppm"%(searchstr,max_val_at/s.get_prop('acq')['SFO1']),human_units=False)
    # }}}
fl.show() # showing the plot
quit()
