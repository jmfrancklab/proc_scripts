""" 
title needed!!
==============

Processes data acquired after applying a triangular pulse to 
a NMR probe.
"""
import numpy as np
from scipy import signal
from pyspecProcScripts import *
from pyspecdata import *
from pyspecProcScripts import postproc_dict
from scipy import signal
import matplotlib.pyplot as plt

fl = figlist_var()
 # {{{ load data, set units, show raw data
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201202_triwave_RMprobe_20ms','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 
    #{{{ raw data
    fl.next('raw together')
    fl.plot(d)
    fl.next('Raw signal for channel 2')
    fl.plot(d['ch',1], alpha=0.5, label='reflection')    
    fl.next('Raw signal for channel 1')
    fl.plot(d['ch',0], alpha=0.5, label='control')
    #}}}
    #{{{pick channel 1 for data and stft signal
    d = d['ch',1]
    d_data = d.data.real
    dw = np.diff(d.getaxis('t')[:2]).item()
    amp=2*sqrt(2)
    fs= 500e3# 1/dw
    f, t, Zxx = signal.stft(d_data, fs, return_onesided=False, nperseg=1000)
    #}}}
    #{{{plotting
    plt.figure()
    x = np.abs(Zxx)
    plt.pcolormesh(t,f,x, vmin=0, vmax=d_data.max(), shading='gouraud')
    plt.title('STFT magnitude for channel 1, 20 ms time/div')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    #}}}
    plt.show()

