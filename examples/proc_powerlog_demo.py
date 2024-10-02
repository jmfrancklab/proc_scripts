import pyspecdata as psd
import pyspecProcScripts as prscr
import os, time, h5py
import matplotlib.pyplot as plt
from sympy import exp as s_exp
from sympy import symbols, Symbol, latex
import numpy as np
from pyspecProcScripts import lookup_table
from Instruments.logobj import logobj
from itertools import cycle

data_target = os.path.normpath(psd.getDATADIR('WK_processed_data'))
filename = "240924_13p5mM_TEMPOL_ODNP_1"
expdata = "ODNP_NMR_comp/ODNP"
nodename = "ODNP"
postproctype = "spincore_ODNP_v5"
my_filename = psd.search_filename(filename+".h5", exp_type=expdata, unique=True)
fl = psd.figlist_var()

Ep_signal_pathway = {'ph1':1}
Ep_f_slice = (-1e3,1e3)
# {{{ this pulls the log from the file
with h5py.File(my_filename,'r') as f:
    log_grp = f['log']
    thislog = logobj.from_group(log_grp)
    log_array = thislog.total_log
    log_dict = thislog.log_dict
# }}}
print(log_array.dtype.names)
# I don't know if the following are code is correct, but you should be able to correct based on the names of the fields, given by the previous line
#log_array['time'] -= log_array['time'][0] # always convert to relative time right away
log_vs_time = psd.nddata(log_array["power"], # new nddata, whose data are the values from the gigatronix
                     [-1], # it's one dimension, whose length is automatically determined
                     ['time'] # the name of hte dimension is time
                     ).setaxis('time', # set the coordinate axis
                               log_array['time'] # to the "time" field of the structured array that comes from the log
                               )
fl.next('power log')
fl.plot(log_vs_time) # should be a picture of the gigatronics powers
# {{{ construct an nddata whose data are the average power values, whose errors are the std of of the power values, and whose time axis is the center time for each power
# {{{ AG does something else, but basically we want to create an nddata that will store our powers and the associated errors
with psd.figlist_var() as fl:
    thisfile, exptype, nodename, post_proc, lookup = (
        filename+".h5",
        expdata,
        nodename,
        postproctype,
        prscr.lookup_table,
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        postproc=post_proc,
        lookup=prscr.lookup_table,
    )
    s.set_units("indirect","s")
#print(s["indirect"])
powername = "time"
dnp_time_axis = s["indirect"]
log_vs_time.set_units("time", 's')
power_vs_time = psd.ndshape([('time',len(dnp_time_axis))]).alloc().set_error(0).set_units('time','s').setaxis('time',np.zeros(len(dnp_time_axis)))
print(dnp_time_axis)
print("log_vs_time is ", log_vs_time)
print(log_vs_time.dimlabels)
relative_times = []
# }}}
print(power_vs_time)
for j,(time_start,time_stop) in enumerate(zip(dnp_time_axis[:]['start_times'],dnp_time_axis[:]['stop_times'])):
    print(log_vs_time['time',-1])
    print(log_vs_time["time":(time_start,time_stop)])
    print(power_vs_time[powername,j])
    power_vs_time[powername,j] = log_vs_time["time":(time_start,time_stop)].mean("time",std=True)
    print(power_vs_time[powername,j])
    fl.next('power log')
    relative_times.append((time_stop-log_vs_time.getaxis('time')[0])+(time_start-log_vs_time.getaxis('time')[0])/2)
    # {{{ these lines show the start and the stop of the power step
    plt.axvline(x=(time_start-log_vs_time.getaxis('time')[0]), color='green', alpha=0.5)
    plt.axvline(x=(time_stop-log_vs_time.getaxis('time')[0]), color='red', alpha=0.5)
    # }}}
print(relative_times)    
power_vs_time.setaxis('time',relative_times)    
fl.plot(power_vs_time, 'o') # this  should be a *single* o at the center of each power step.  Its y value should be the avaerage power for that step, and its error bars should give the standard deviation of the power over the step
# }}}
fl.show()
