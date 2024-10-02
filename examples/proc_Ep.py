"""
Process Enhancement experiment 
====================================================
Opens .h5 results file, uses rough_table_of_integrals() to roughly process dataset including generating a table of integrals
"""


import pyspecProcScripts as prscr
import pyspecdata as psd
import numpy as np
import os, time, h5py
import matplotlib.pyplot as plt
from sympy import exp as s_exp
from sympy import symbols, Symbol, latex
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

with psd.figlist_var() as fl:
    thisfile, exptype, nodename, post_proc, lookup = (
        "240924_13p5mM_TEMPOL_ODNP_1.h5",
        "ODNP_NMR_comp/ODNP",
        "ODNP",
        "spincore_ODNP_v5",
        prscr.lookup_table,
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        postproc=post_proc,
        lookup=prscr.lookup_table,
    )
    s["indirect"] = s["indirect"]["start_times"]
    s.set_units("indirect", "s")
    prscr.rough_table_of_integrals(s, fl=fl)
    
    # {{{ Power log [probably want as separate function?]
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
    # {{{ Initialize needed variables
    Ep_signal_pathway = {'ph1':1}
    Ep_f_slice = (-1e3,1e3)
    powername = "indirect"
    # }}}
    # {{{ Offset, phasing, and alignment --- it seems like some/all of this is unnecessary    
    s.ift('t2')
    s.set_units('t2','s')
    s /= prscr.zeroth_order_ph(prscr.select_pathway(s['t2':0],Ep_signal_pathway))
    s_align = s.C
    FID = s.C
    FID = FID['t2':(0,None)]
    FID *= 2
    FID['t2':0] *= 0.5
    FID.ft('t2')
    s_align = s_align['t2':(0,None)]
    s_align *= 2
    s_align['t2':0] *= 0.5
    s_align.ft('t2')
    s_align.ift('t2')
    fl.next('test time')
    fl.text('I am making a copy of the data to apply apodization and then obtain the optimal shift from the alignment. But this apodization is not applied to the actual data.')
    filter_timeconst = 32e-3
    myfilter = np.exp(-abs((s_align.fromaxis('t2')-s_align.get_prop('acq_params')['tau_us']*1e-6))/filter_timeconst)
    fl.plot(myfilter * abs(s_align.C).max())
    s_align *= myfilter
    s_align.ft('t2')
    if 'nScans' in s.dimlabels:
        mysgn = prscr.determine_sign(prscr.select_pathway(s_align['t2':(Ep_f_slice[0]+500,Ep_f_slice[-1]-500)].C.mean('nScans'),Ep_signal_pathway), (Ep_f_slice[0]+500,Ep_f_slice[-1]-500))
    else:
        mysgn = prscr.determine_sign(prscr.select_pathway(s_align['t2':(Ep_f_slice[0]+500,Ep_f_slice[-1]-500)],Ep_signal_pathway), (Ep_f_slice[0]+500,Ep_f_slice[-1]-500))
    s_align *= mysgn
    s_align.reorder(['t2'],first = False)
    s_align.ift(['ph1'])
    s_align.reorder(['t2'],first = False)
    this_sigma = 10
    s_align.set_error(None)
    opt_shift,sigma, my_mask = prscr.correl_align((s_align.C),
            indirect_dim=powername,
            signal_pathway=Ep_signal_pathway,
            sigma = this_sigma, fl=fl)
    s.ft('t2')
    s.ift(['ph1'])
    s.ift('t2')
    s *= np.exp(-1j*2*np.pi*opt_shift*s.fromaxis('t2'))
    s.ft(['ph1'])
    # }}}
    # {{{ FID slice
    d=s.C
    d = d['t2':(0,None)]
    d *= 2
    d['t2':0] *= 0.5
    d.ft('t2')
    print(d)
    fl.next('E(p) FID sliced and time sliced after alignment')
    if 'nScans' in s.dimlabels:
        fl.image(d.C.mean('nScans'))
    else:
        fl.image(d)
    # }}}
    #{{{Integrate with error
    #{{{finding average power over steps
    dnp_time_axis = d.C.getaxis(powername).copy()
    print("dnp_time_axis is ", dnp_time_axis)
    print(log_start_time)
    print(dnp_time_axis[:]['start_times'])
    dnp_time_axis[:]['start_times'] -= log_start_time
    dnp_time_axis[:]['stop_times'] -= log_start_time
    nddata_time_axis = nddata(dnp_time_axis,[-1],[powername])
    new_time_axis = nddata_time_axis.C.data
    new_time_axis = nddata(new_time_axis,[-1],[powername])
    power_vs_time = ndshape(nddata_time_axis).alloc().set_units(powername,'s')
    power_vs_time.set_error(0)
    power_vs_time.setaxis(powername,new_time_axis.data)
    #}}}
    #{{{find values for Ep
    for j,(time_start,time_stop) in enumerate(zip(dnp_time_axis[:]['start_times'],dnp_time_axis[:]['stop_times'])):
        power_vs_time[powername,j] = power_axis[powername:((time_start),(time_stop))].mean(powername,std=True)
        power_vs_time.set_units(powername,'s')
        fl.next('Instrument Power Log')
        plt.axvline(x=time_start,color='red',alpha=0.5)
        plt.axvline(x=time_stop,linestyle=':',color='red',alpha=0.5)
    avg_p_vs_t = nddata(power_vs_time.data,[-1],[powername])
    avg_p_vs_t.set_error(power_vs_time.get_error())
    avg_p_vs_t.setaxis(powername,dnp_time_axis[:]['start_times'])
    avg_p_vs_t.data[0] = 0
    avg_p_vs_t[powername,0].set_error(0)
    fl.plot(avg_p_vs_t,'ro',capsize=6)
    plt.xlabel('Time / s')
    plt.ylabel('Power / W')
    power_axis_Ep = np.real(power_vs_time.data)
    #}}}
    #d.mean('nScans')
    d.rename(powername,'power')
    d.setaxis('power',power_axis_Ep)
    fl.next('Thermals before integration')
    for j in range(thermal_scans):
        fl.plot(select_pathway(d['power',j].C,Ep_signal_pathway),label='thermal %d'%j)
    mysgn = determine_sign(select_pathway(d['t2':(Ep_f_slice[0]+500,Ep_f_slice[-1]-500)].C,Ep_signal_pathway), (Ep_f_slice[0]+500,Ep_f_slice[-1]-500))
    d *= mysgn
    error_pathway = (set(((j) for j in range(ndshape(d)['ph1'])))
            - set(excluded_pathways)
            -set([(Ep_signal_pathway['ph1'])]))
    error_pathway = [{'ph1':j} for j in error_pathway]
    s_int,this_frq_slice = integral_w_errors(d.C,Ep_signal_pathway,error_pathway,
            indirect='power',cutoff = 0.01,return_frq_slice=True)
    fl.next('E(p) Integration Limits')
    fl.plot(select_pathway(d,Ep_signal_pathway))
    plt.axvline(x=this_frq_slice[0])
    plt.axvline(x=this_frq_slice[1])
    #}}}
    
