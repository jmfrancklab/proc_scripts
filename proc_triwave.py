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
    fl.next('Raw signal %s'%searchstr)
    fl.plot(d['ch',0], alpha=0.5, label='control')    
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    #fl.show();quit()
    d_data = d.data
    dw = np.diff(d.getaxis('t')[:2]).item()
    amp=2*sqrt(2)
    fs = 50e6
    f, t, Zxx = signal.stft(d_data, fs, nperseg=100)
    np.reshape(Zxx,(3,2))
    print("shapes of t, f, Z", t.shape,f.shape,Zxx.shape)
    #quit()
    x = np.abs(Zxx)
    plt.pcolormesh(t,f,x,cmap='RdYlGn')
    plt.title('STFT magnitude')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.show()

