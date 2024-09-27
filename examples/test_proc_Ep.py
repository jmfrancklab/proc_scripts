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

data_target = os.path.normpath(psd.getDATADIR('WK_processed_data'))
fl = psd.figlist_var()
Ep_signal_pathway = {'ph1':1}
Ep_f_slice = (-1e3,1e3) # stop using this, switch to peakrange (?)
filename = '240924_13p5mM_TEMPOL_ODNP_1'
exptype = 'ODNP_NMR_comp/ODNP'
postproctype = 'spincore_ODNP_v3'
nodename = 'ODNP'
fl.basename = nodename
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
my_filename = psd.search_filename(filename+".h5", exp_type=exptype, unique=True)

# {{{ power log
powername = 'time'
with h5py.File(my_filename,'r') as f:
    log_grp = f['log']
    thislog = logobj.from_group(log_grp)
    read_array = thislog.total_log
    read_dict = thislog.log_dict
for j in range(len(read_array)):
    thistime,thisrx,thispower,thiscmd = read_array[j]
log_start_time = read_array[powername][0]
relative_time = read_array[powername]
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
    ax_Rx.axvline(x=thisevent[powername]-log_start_time)
    ax_power.axvline(x=thisevent[powername]-log_start_time)
    y_pos = j/n_events
plt.tight_layout()
power_axis = psd.nddata(read_array['power'],[-1],[powername])
power_axis.setaxis(powername,relative_time)
power_axis.setaxis(powername,lambda x: x - log_start_time)
power_axis = psd.nddata(read_array['power'],[-1],[powername])
power_axis.setaxis(powername,relative_time)
power_axis.name('power')
power_axis.data = (10**((power_axis.data+coupler_atten)/10+3))/1e6 #convert to W
fl.next('Instrument Power Log')
fl.plot(power_axis,'.')
#fl.show();quit()
# }}}
# {{{ raw proc E(p)
s = psd.find_file(filename, exptype, expno=nodename, postproc=postproctype, lookup=lookup_table)
fl.next('Raw E(p)',figsize = (5,20))
fl.image(s)
print("s is ", s)
print("shape of s is ", psd.ndshape(s))
# }}} 
# {{{ DC offset correction 
s.ift('t2')
s.ift(['ph1'])
Ep_t_max = s.getaxis('t2')[-1]
rx_offset_corr = s['t2':(Ep_t_max*0.75, None)]
rx_offset_corr = rx_offset_corr.mean(['t2'])
s -= rx_offset_corr
s.ft('t2')
s.ft(['ph1'])
# }}} 
s = prscr.rough_table_of_integrals(s, fl=fl, title=conctag)
fl.show()
