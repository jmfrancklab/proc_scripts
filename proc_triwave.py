import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from proc_scripts import *
from pyspecdata import *
from proc_scripts import postproc_dict


fl = figlist_var()
 # {{{ load data, set units, show raw data
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201103_triwave_coile_1','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 
    fl.next('Raw signal %s'%searchstr)
    fl.plot(d['ch',0], alpha=0.5, label='control')    
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    fl.show();quit()
    T = 10 #time duration in seconds
    fs = 1000020 # sampling rate in Hz
    fm = 20000 #frequency of modulation in Hz
    B = 1000000 #bandwidth in Hz
    frame_size = 0.000010 #frame size in seconds
    hop = 0.000005 #hop size

    t = np.linspace(0,T,T*fs,endpoint=False)
    dt = t[1]-t[0]
    f = simu_freq_sawtooth(t, fm, B, fd=0, width=0.5)
    y = d
    Y = stft(y,fs,frame_size,hop)

    t_min = t[0]
    t_max = t[-1]
    f_min = 0
    f_max = int(fs * 0.3)
    plt.figure(figsize=(20,15))

    ax = plt.subplot(2,1,1)
    plt.plot(t,f)
    plt.title('signal')
    plt.ylim(t_min,t_max)
    plt.ylim(f_min,f_max)
    plt.ylabel('frequency(Hz)')
    ax = plt.subplot(2,1,2)
    Fz = int(frame_size * fs * 0.3)
    ax.imshow(np.absolute(Y[:, :Fz].T), origin='lower',
            aspect='auto', interpolation='nearest',
            extent=[t_min,t_max,f_min,f_max])
    plt.title('measured signal')
    plt.xlabel('time(s)')
    plt.ylabe('frequency in Hz')
    fl.show();quit()
