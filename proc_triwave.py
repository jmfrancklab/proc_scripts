import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from proc_scripts import *
from pyspecdata import *
from proc_scripts import postproc_dict


fl = figlist_var()
 # {{{ load data, set units, show raw data
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201113_triwave_control_2','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 
    fl.next('Raw signal %s'%searchstr)
    fl.plot(d['ch',0], alpha=0.5, label='control')    
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    fl.show();quit()
    d.stft(d,50e6,0.5,0.25)
    fl.next('after STFT')
    fl.plot(d['ch',0],alpha=0.5,label='control')
    fl.plot(d['ch',1],alpha=0.5,label='reflection')
    fl.show();quit()

