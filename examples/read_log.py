from pyspecdata import *
from pyspecProcScripts import logobj
import h5py
from matplotlib.pyplot import subplots
from numpy import log10
coupler_atten = 22
myfilename = search_filename('230626_batch230515_E37_Ras_B10_ODNP_1.h5',exp_type='ODNP_NMR_comp/ODNP',unique = True)
with figlist_var() as fl:
    # {{{ open h5 file to real log
    with h5py.File(myfilename,'r') as f:
        log_grp = f['log']
        thislog = logobj.from_group(log_grp)
        read_array = thislog.total_log
        read_dict = thislog.log_dict
    # }}}
    # {{{ separate the time, RX, power, and commands from the 
    #     log array
    for j in range(len(read_array)):
        thistime,thisrx,thispower,thiscmd = read_array[j]
    # }}}
    # In order to properly set the time axis to start at 0
    # both the log's start time will be subtracted from the
    # the relative time recorded
    log_start_time = read_array['time'][0]
    relative_time = read_array['time']
    time_axis = relative_time - log_start_time
    # }}}
    # {{{ plot the output power and reflection
    fig, (ax_Rx,ax_power) = subplots(2,1,figsize=(10,8))
    fl.next('log figure', fig = fig)
    ax_Rx.set_ylabel('Rx / mV')
    ax_Rx.set_xlabel('Time / ms')
    ax_Rx.plot(time_axis, read_array['Rx'],'.')
    ax_power.set_ylabel('power / W')
    ax_power.set_xlabel('Time / ms')
    ax_power.plot(time_axis,10**(read_array['power']/10 + 3+coupler_atten/10),'.')
    # }}}
    # {{{ Add a vertical line at the time the data acquisition for the
    #     set power began
    mask = read_array['cmd'] != 0
    n_events = len(relative_time[mask])
    for j,thisevent in enumerate(read_array[mask]):
        ax_Rx.axvline(x=thisevent['time'] - log_start_time)
        ax_power.axvline(x=thisevent['time']-log_start_time)
    # }}}
    plt.tight_layout()
