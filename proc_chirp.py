from pyspecdata import *
from scipy.optimize import minimize,leastsq
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
from scipy.signal import tukey
do_slice = False # slice frequencies and downsample -- in my hands, seems to decrease the quality of the fit 
standard_cost = False # use the abs real to determine t=0 for the blip -- this actually doesn't seem to work, so just use the max
show_transfer_func = False # show the transfer function -- will be especially useful for processing non-square shapes
logger = init_logging('info')
#init_logging(level='debug')
# 2to3 JF 1/31

fl = figlist_var()
 # {{{ load data, set units, show raw data
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201026_chirp_cap_probe_2','test_equip','capture1','chirp',True)
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 
    print(d.getaxis('t'))
    fl.next('Raw signal %s'%searchstr)
    fl.plot(d['ch',0], alpha=0.5, label='control')    
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    d.ft('t',shift=True)
    fl.next('freq domain for coil sphere removed untuned')
    fl.plot(d['ch',0],alpha=0.5,label='control')
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    d_control = d['ch',0]['t':(0,None)]
    d_refl = d['ch',1]['t':(0,None)]
    d = d_refl/d_control
    fl.next('c')
    fl.plot(abs(d))
    ylim(0,1)
    fl.show();quit()

