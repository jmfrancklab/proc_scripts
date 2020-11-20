import numpy as np
from scipy import signal
from proc_scripts import *
from pyspecdata import *
from proc_scripts import postproc_dict
from scipy import signal
import matplotlib.pyplot as plt



fl = figlist_var()
 # {{{ load data, set units, show raw data
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201113_triwave_control_2','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 
    fl.next('Raw signal for channel 1')
    fl.plot(d['ch',1], alpha=0.5, label='control')    
    fl.next('Raw signal for channel 0')
    fl.plot(d['ch',0], alpha=0.5, label='reflection')
    d = d['ch',1]
    d_data = d.data.real
    print("shape of d before stft",d_data.shape)
    dw = np.diff(d.getaxis('t')[:2]).item()
    amp=2*sqrt(2)
    fs= 1/dw
    f, t, Zxx = signal.stft(d_data, fs, return_onesided=True, nperseg=1000)
    print("shapes of t, f, Z", t.shape,f.shape,Zxx.shape)
    plt.figure()
    x = np.abs(Zxx)
    plt.pcolormesh(t,f,x, vmin=0, vmax=d_data.max(), shading='gouraud')
    plt.title('STFT magnitude')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.show()

