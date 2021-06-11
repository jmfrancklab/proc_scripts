"""
Short Time Fourier Transform functions
======================================

taken from 
http://tsaith.github.io/time-frequency-analysis-with-short-time-fourier-transform.html

"""
#http://tsaith.github.io/time-frequency-analysis-with-short-time-fourier-transform.html
from pyspecdata import *
from pyspecProcScripts import *
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (20,15)
plt.rcParams['image.interpolation'] = 'nearest'
plt.rcParams.update({'font.size': 22})

def stft(x, fs, frame_size, hop):
    """
    Perform STFT (short-time fourier transform).

    x:Input Data.
    fs: sampling rate.
    frame_size: frame size.
    hop: hop size
    """
    frame_samp = int(frame_size*fs)
    hop_samp = int(hop*fs)
    w = np.hanning(frame_samp) #Hanning window
    X = np.array([np.fft.fft(w*x[i:i+frame_samp])
        for i in range(0, len(x)-frame_samp,
            hop_samp)])
    return X

def isftf(X,fs,T,hop):
    """
    Perform inverse STFT.

    X: Input data.
    fs: sampling rate.
    T: total time duration.
    hop: Hop size.
    """
    x = np.zeros(T*fs)
    frame_samp = X.shape[1]
    hop_samp = int(hop*fs)
    for n,i in enumerate(range(0, len(x)-frame_sampe,
        hop_samp)):
        x[i:i+frame_samp] += np.real(np.fft.ifft(X[n]))
    return x

def simu_waves(f,dt,amp0=1,phi0=0):
    """
    Return the simulated waves.
    y(t) = amp0 * cos(phi(t) + phi0),
    where phi(t) = 2*pi*\int_0^t f(t) * dt.

    f: Instantaneous frequencies.
    dt: time interval
    amp0: amplitude
    phi0: initial phase. when it is -pi/2, sin waves 
    are produced.
    """
    phi = 2 * np.pi * np.cumsum(f) * dt
    y = amp0*np.cos(phi + phi0)
    return y

def simu_freq_sawtooth(t, fm=1, B=1, fd=0, width=0.5):
    '''
    simulated frequencies of sawtooth modulation.

    t: time array
    fm: modulation frequency
    fd: doppler frequency shift
    B: bandwidth
    '''
    f = B*0.5*(signal.sawtooth(2 * np.pi * fm * t,
        width=width) +1)
    f += fd
    return f
