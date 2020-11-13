import numpy as np
import matplotlib.pyplot as plt
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
    print(nddata(d_data))
    dw = diff(d.getaxis('t')[:2]).item()
    amp=2*sqrt(2)
    fs = 50e6
    f, t, Zxx = signal.stft(d_data, fs, nperseg=1000)
    plt.pcolormesh(t, f, np.abs(Zxx), vmin=0, vmax=amp, shading='gouraud')
    plt.title('STFT magnitude')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.show()

