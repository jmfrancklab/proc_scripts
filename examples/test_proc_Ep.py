import pyspecdata as psd
import pyspecProcScripts as prscr
import numpy as np
import matplotlib.pyplot as plt
import os, time, h5py
from sympy import exp as s_exp
from sympy import symbols, Symbol, latex
from numpy import exp
from pyspecProcScripts import lookup_table
from Instruments.logobj import logobj
from itertools import cycle

data_target = os.path.normpath(getDATADIR('WK_processed_data'))
fl = figlist_var()
Ep_signal_pathway = {'ph1':1}
Ep_f_slice = (-1e3,1e3) # stop using this, switch to peakrange (?)
thermal_scans = 2 # pull this from the h5 file  
filename = '240924_13p5mM_TEMPOL_ODNP_1'
exptype = 'ODNP_NMR_comp/ODNP'
postproctype = 'spincore_ODNP_v3'
powers = []
errors = []
power_list = []
start_times = []
coupler_atten = 22 # can I pull this from h5?
excluded_pathways = [(0,0),(0,3)] 
filename_out = 'VTprobe_test_240924.h5' 
conctag = '240924_13p5mM_TEMPOL_ODNP_1' 
color_cycle = cycle(['red','orange','yellow','green','cyan','blue','purple','magenta',
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
    '#7f7f7f', '#bcbd22', '#17becf'])
my_filename = search_filename(filename+".h5", exp_type=exptype, unique=True)
# {{{ power log
with h5py.File(my_filename,'r') as f:
    log_grp = f['log']
    thislog = logobj.from_group(log_grp)
    read_array = thislog.total_log
    read_dict = thislog.log_dict
for j in range(len(read_array)):
    thistime,thisrx,thispower,thiscmd = read_array[j]
log_start_time = read_array['time'][0]
relative_time = read_array['time']
fig, (ax_Rx,ax_power) = plt.subplots(2,1, figsize=(10,8))
fl.next("log figure",fig=fig)
ax_Rx.set_ylabel('Rx / mV')
ax_Rx.set_xlabel('Time / ms')
ax_Rx.plot(relative_time-log_start_time, read_array['Rx'],'.')
ax_power.set_ylabel('power / W')
ax_power.set_xlabel('Time / ms')
ax_power.plot(relative_time-log_start_time,10**(read_array['power']/10+3+coupler_atten/10),'.')
mask = read_array['cmd'] != 0
n_events = len(relative_time[mask])
for j, thisevent in enumerate(read_array[mask]):
    ax_Rx.axvline(x=thisevent['time']-log_start_time)
    ax_power.axvline(x=thisevent['time']-log_start_time)
    y_pos = j/n_events
plt.tight_layout()
power_axis = nddata(read_array['power'],[-1],['time'])
power_axis.setaxis('time',relative_time)
power_axis.setaxis('time',lambda x: x - log_start_time)
power_axis = nddata(read_array['power'],[-1],['time'])
power_axis.setaxis('time',relative_time)
power_axis.name('power')
power_axis.data = (10**((power_axis.data+coupler_atten)/10+3))/1e6 #convert to W
fl.next('Instrument Power Log')
fl.plot(power_axis,'.')
powername = 'time'
#fl.show();quit()
# }}}
# {{{ proc E(p)
s = find_file(filename, exptype, expno='ODNP', postproc=postproc
              postproc=postproctype, lookup=lookup_table)
if 'indirect' in s.dimlabels:
        s.rename('indirect','time')
    if 'power' in s.dimlabels:
        s.rename('power','time')
    if 'nScans' in s.dimlabels: 
        s.reorder(['nScans','ph1',powername])
    else:
        s.reorder(['ph1',powername])
    fl.next('Raw E(p)',figsize = (5,20))
    fl.image(s)
fl.show();quit()
# }}}
