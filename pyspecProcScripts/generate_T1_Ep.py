import os, time, h5py
import pylab as plt
from numpy import empty,exp
from sympy import exp as s_exp
from sympy import symbols, Symbol, latex
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from Instruments.logobj import logobj
from itertools import cycle
fl = figlist_var()
color_cycle = cycle(['red','orange','yellow','green','cyan','blue','purple','magenta',
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
    '#7f7f7f', '#bcbd22', '#17becf'])
def generate_T1_Ep(filename,
        has_Ep = True,
        has_IR = True,
        Ep_signal_pathway={'ph1':1},
        IR_signal_pathway = {'ph1':0,'ph2':-1},
        Ep_nodename = 'ODNP',
        IR_nodenames = None,
        IR_postproc='spincore_IR_v1',
        Ep_postproc = 'spincore_ODNP_v4',
        Ep_f_slice = None,
        IR_f_slice = None,
        Ep_apo = 50e-3,
        IR_f_slice_auto_scale = 500,
        Ep_f_slice_auto_scale = 500,
        Mi_max = 2e6,
        drift_max = 200,
        Ep_drift_max = 200,
        nPowers = 15,
        W=None,
        IR_cutoff = 0.15,
        Ep_cutoff = 0.15,
        extra_IR_zero = False,
        Ep_extra_zero = False,
        log=True,
        Ep_flip = False,
        fl=None):
    T1_list = []
    errors = []
    power_list = []
    start_times = []
    powers = []
    excluded_pathways = [(0,0)]
    coupler_atten = 22
    if log:
        #{{{load in log
        my_filename = search_filename(filename+".h5",exp_type='ODNP_NMR_comp/ODNP',unique=True)
        with h5py.File(my_filename,'r') as f:
            log_grp = f['log']
            thislog = logobj.from_group(log_grp)
            read_array = thislog.total_log
            read_dict = thislog.log_dict
        for j in range(len(read_array)):
            thistime,thisrx,thispower,thiscmd = read_array[j]
        log_start_time = read_array['time'][0]
        relative_time = read_array['time']
        if fl is not None:
            fig, (ax_Rx,ax_power) = plt.subplots(2,1, figsize=(10,8))
            fl.next("log figure",fig=fig)
            ax_Rx.set_ylabel('Rx/mV')
            ax_Rx.plot(relative_time-log_start_time, read_array['Rx'],'.')
            ax_power.set_ylabel('power/W')
            ax_power.plot(relative_time-log_start_time,10**(read_array['power']/10+3+coupler_atten/10),'.')
        mask = read_array['cmd'] != 0
        n_events = len(relative_time[mask])
        if fl is not None:
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
        if fl is not None:
            fl.next('Instrument Power Log')
            fl.plot(power_axis,'.')
        #}}}
        powername = 'time'
    else:
        powername = 'time'
    if has_IR:
        if fl is not None:
            # have these come right after the log
            fl.text("Overall fits:\n\n")
            fl.next('IR fit - before norm')
            fl.next('IR fit - normalized')
            fl.text(r"\par\par")
    #{{{IR processing
        if IR_nodenames is None:
            # if set to None, just assume any node with
            # "IR" is an IR -- this could be misleading!
            #
            # note that this requires post 7/11/22
            # pyspecdata
            all_nodes = find_file(filename,
                    exp_type='ODNP_NMR_comp/ODNP',
                    return_list=True)
            IR_nodenames = [j for j in all_nodes if
                    'IR' in j]
        for (j,nodename) in enumerate(IR_nodenames):
            fl.text("processing "+lsafen(nodename))
            IR = find_file(filename, exp_type='ODNP_NMR_comp/ODNP', expno=nodename,
                    postproc=IR_postproc, lookup=lookup_table)
            #{{{using manual powers
            if manual_IR_meter_powers is not None:
                try:
                    IR.get_prop('acq_params')['meter_power'] = manual_IR_meter_powers[nodename]
                except:
                    raise ValueError(f"you need to supply the power for node {nodename}")
            if hack_IR_powers:
                if('meter_power' in IR.get_prop('acq_params').keys()): raise ValueError("you don't need to hack your powers!")
                m = hack_IR_nopower_re.match(nodename)
                if m is not None:
                    IR.get_prop('acq_params')['meter_power'] = 0
                else:
                    m = hack_IR_re.match(nodename)
                    if m is None: raise ValueError(f"your node name {nodename} doesn't match the usual pattern for hacking!")
                    IR.get_prop('acq_params')['meter_power'] = 10**(0.1*np.double(m.groups()[0])-3)
            #}}}
            if older:
                IR = IR['indirect',j]
            else:
                pass
            IR.reorder(['ph1','ph2','nScans','vd'])
            IR['ph2',0]['ph1',0]['t2':0]=0 #kill axial noise
            #{{{DC offset correction
            IR.ift('t2')
            IR.ift(['ph1','ph2'])
            t_max = IR.getaxis('t2')[-1]
            rx_offset_corr = IR['t2':(t_max*0.75,None)]
            rx_offset_corr = rx_offset_corr.mean(['t2'])
            IR -= rx_offset_corr
            IR.ft('t2')
            IR.ft(['ph1','ph2'])
            #}}}
            if IR_f_slice is None:
                for_lims = IR.C.mean('nScans')
                drift_bounds = [
                        abs(select_pathway(for_lims['vd',0],IR_signal_pathway)).argmax().item()
                        ,
                        abs(select_pathway(for_lims['vd',-1],IR_signal_pathway)).argmax().item()
                        ]
                drift_bounds.sort()
                this_IR_f_slice = array(drift_bounds) + r_[-IR_f_slice_auto_scale,IR_f_slice_auto_scale]
            else:
                this_IR_f_slice = IR_f_slice
            zero_crossing = abs(select_pathway(IR['t2':this_IR_f_slice].C.mean('nScans'),IR_signal_pathway)).C.sum('t2').argmin('vd',raw_index=True).item()
            if fl is not None:
                fl.next('Raw IR for %s'%nodename)
                fl.image(IR.C.mean('nScans'))
                #fl.show();quit() #WARREN,UNCOMMENT THIS THE FIRST TIME YOU RUN SO YOU CAN SEE WHAT YOUR IR_f_slice VALUES SHOULD BE
            IR=IR['t2':this_IR_f_slice]
            IR.ift('t2')
            ##{{{phasing
            IR.set_units('t2','s')
            IR.ft('t2')
            IR.ift('t2')
            IR.set_units('t2','s')
            best_shift,cost_fn = hermitian_function_test(select_pathway(IR.C.mean('vd'),
                IR_signal_pathway),
                echo_before=IR.get_prop('acq_params')['tau_us']*1e-6*1.5,
                basename=f"Herm for {nodename}",fl=fl)
            print("THE ESTIMATED TAU IS %.6f"%best_shift)
            actual_tau = IR.get_prop('acq_params')['tau_us']/1e6
            if (best_shift < actual_tau-0.8e-3) or (best_shift > actual_tau+0.8e-3):
                if fl is not None:
                    fl.text(r'\textcolor{red}{\textbf{I am hard-setting the first-order phase for dataset %s}}'%nodename)
                    print("HARD SETTING")
                    fl.basename = nodename
                    best_shift = hermitian_function_test(select_pathway(IR.C,
                        IR_signal_pathway),
                        echo_before=IR.get_prop('acq_params')['tau_us']*1e-6*1.5,
                        )
                    fl.basename=None
                best_shift = actual_tau + 0.0002
            IR.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
            IR /= zeroth_order_ph(select_pathway(IR['t2',0],IR_signal_pathway))
            IR.ft('t2')
            if fl is not None:
                fl.next('phased and FID sliced')
                fl.image(IR.C.mean('nScans'),human_units=False)
                #fl.show();quit()
             #}}}
            #{{{Alignment
            first = abs(select_pathway(IR['vd',0].C.mean('nScans'),IR_signal_pathway)).argmax().item()
            last =abs(select_pathway(IR['vd',-1].C.mean('nScans'),IR_signal_pathway)).argmax().item() 
            if last > first:
                drift = last - first
            else:
                drift = first - last
            print("Your signal has a drift of %.4f"%drift)    
            if drift < drift_max:
                IR.ift(['ph1','ph2'])
                IR['vd',:zero_crossing] *= -1
                opt_shift,sigma, my_mask = correl_align(IR.C,
                        indirect_dim='vd',
                        signal_pathway=IR_signal_pathway,
                        sigma=1500,
                        )
                IR.ift('t2')
                IR *= np.exp(-1j*2*pi*opt_shift*IR.fromaxis('t2'))
                IR.ft(['ph1','ph2'])
                IR.ft('t2')
                IR['vd',:zero_crossing] *= -1
            else:
                drift /= len(IR.getaxis('vd'))
                IR.ift('t2')
                for i in range(len(IR.getaxis('vd'))):
                    IR['vd',i] *= np.exp(-1j*2*pi*(drift+i*drift)*IR.fromaxis('t2'))
                IR.ft('t2')
                IR['vd',:zero_crossing] *= -1
                ph0 = IR['t2':(-900,-500)].C.sum('t2')
                ph0 /= abs(ph0)
                IR /= ph0
                IR.ift(['ph1','ph2'])
                opt_shift,sigma, my_mask = correl_align((IR),indirect_dim='vd',
                        signal_pathway=IR_signal_pathway,sigma = 1500)
                IR.ift('t2')
                IR *= np.exp(-1j*2*pi*opt_shift*IR.fromaxis('t2'))
                IR.ft('t2')
                IR.ft(['ph1','ph2'])
            #}}}
            IR.ift('t2')
            #{{{FID slice
            d=IR.C
            d = d['t2':(0,None)]
            d['t2':0] *= 0.5
            d.ft('t2')
            if fl is not None:
                fl.next('IR aligned and FID sliced - %s'%nodename)
                fl.image(d.C.mean('nScans'),human_units=False)
                #fl.show();quit() #WARREN, UNCOMMENT THE HASH TO SEE THE ALIGNED DATA,IF THIS DOESN'T WORK, TRY CHANGING THE KWARG "drift_max" TO BE HIGHER OR LOWER THAN THE SIGNAL DRIFT PRINTED IN THE TERMINAL
            #}}}
            #{{{Integrate with error
            if d['vd',0] < 0:
                d['vd',:zero_crossing] *= -1
            this_IR = (d.C.mean('nScans'))
            if extra_IR_zero:
                d.ift('t2')
                d /= zeroth_order_ph(select_pathway(d['t2',0],IR_signal_pathway))
                d.ft('t2')
            error_pathway = (set(((j,k) for j in range(ndshape(d)['ph1']) for k in range(ndshape(d)['ph2'])))
                    - set(excluded_pathways)
                    -set([(IR_signal_pathway['ph1'],IR_signal_pathway['ph2'])]))
            error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
            s_int,frq_slice = integral_w_errors(abs(this_IR.C),IR_signal_pathway,error_pathway,
                    cutoff = IR_cutoff,
                    indirect='vd',return_frq_slice=True)
            if fl is not None:
                fl.next('IR Integration Limits - %s'%nodename)
                for j in range(len(d.getaxis('vd'))):
                    fl.plot(select_pathway(this_IR['vd',j],
                        IR_signal_pathway),label = '%d'%j,human_units=False)
                plt.axvline(frq_slice[0])
                plt.axvline(frq_slice[-1])
                #fl.show();quit()#WARREN, UNCOMMENT THIS TO SEE IF THE INTEGRATION LIMITS LOOK OKAY, YOU CAN ADJUST THE KWARG IR_cutoff TO MOVE THESE AROUND
            #}}}
            #{{{Fitting
            thiscolor = next(color_cycle)
            if (s_int['vd',0]>0) and (s_int['vd',1]>0):
                print("HEY")
                s_int['vd',:zero_crossing] *= -1
            elif (s_int['vd',0]>0):
                print("HEY HEY")
                s_int['vd',:zero_crossing] *= -1
            M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd",real=True)
            logger.debug(strm("acq keys",IR.get_prop('acq_params')))
            if W is None:
                W = (IR.get_prop('acq_params')['FIR_rep']*1e-6
                        +
                        IR.get_prop('acq_params')['acq_time_ms']*1e-3)
            else:
                W = W
            functional_form = Mi*(1-(2-s_exp(-W*R1))*s_exp(-vd*R1))
            IR_data = nddata(s_int.data,['vd'])
            IR_data.setaxis('vd',s_int.getaxis('vd'))
            f = lmfitdata(IR_data)
            f.functional_form = functional_form
            f.set_guess(
                    M_0 = dict(value=-Mi_max+1, min = -Mi_max, max = 0),
                    M_inf = dict(value=Mi_max-1, min = 0, max=Mi_max),
                    R_1 = dict(value=5, min=0.001, max=100)
                    )
            f.settoguess()
            guess = f.eval(100)
            f.fit()
            T1 = 1./f.output('R_1')
            Mi = f.output('M_inf')
            T1_list.append(T1)
            this_ls = "-"
            fit = f.eval(100)
            print(T1)
            if fl is not None:
                fl.next('IR fit - before norm')
                fl.plot(s_int, 'o', color=thiscolor, label='%s'%nodename)
                fl.plot(fit, ls = this_ls, color=thiscolor,alpha=0.5,
                        label = 'fit for %s'%nodename)
                fl.next('IR fit - normalized')
                fl.plot(s_int/Mi, 'o', color=thiscolor, label=nodename)
                fl.plot(fit/Mi, 
                        ls=this_ls,
                        color=thiscolor,
                        alpha=0.5,
                        label='fit for %s'%nodename,
                        )
                ax=plt.gca()
                #fl.show();quit()
            #}}}
            #{{{set power axis
            if log:
                start_time = IR.get_prop('start_time')-log_start_time
                stop_time = IR.get_prop('stop_time')-log_start_time
                avg_power = power_axis['time':(start_time,
                    stop_time)].mean('time',std=True)
                errors.append(avg_power.get_error())
                power_list.append(avg_power)
                powers.append(avg_power.data)
                start_times.append(start_time)
                if fl is not None:
                    fl.next('Instrument Power Log')
                    plt.axvline(x=start_time,color='k',alpha=0.5)
                    plt.axvline(x=stop_time,linestyle = ':', color='k',alpha=0.5)
                nddata_p_vs_t = nddata(power_list,[-1],['time'])
                nddata_p_vs_t.setaxis('time',start_times)
                for j in range(len(power_list)):
                    nddata_p_vs_t['time',j] = power_list[j]
                nddata_p_vs_t.set_error(errors)
                if fl is not None:
                    fl.plot(nddata_p_vs_t,'ko',capsize=6)
            #}}}        
        #{{{make T1p 
        T10 = T1_list[0]
        T1p = nddata(T1_list,[-1],['power'])
        T1p.setaxis('power',powers)
        if fl is not None:
            fl.next(r'$T_{1}$(p)')
            fl.plot(T1p,'o')
            plt.ylabel(r'$T_{1}$ / s')
            plt.xlabel('Power / W')
            #fl.show();quit()#WARREN UNCOMMENTING THIS WILL LET YOU SEE THE FIT! 
    #}}}
    #}}}
    if has_Ep:
        #{{{Ep processing    
        s = find_file(filename, exp_type = 'ODNP_NMR_comp/ODNP', expno= Ep_nodename,
                postproc=Ep_postproc, lookup = lookup_table)
        if 'indirect' in s.dimlabels:
            s.rename('indirect','time')
        if 'power' in s.dimlabels:
            s.rename('power','time')
        s.reorder(['ph1',powername])
        if 'nScans' in s.dimlabels:
            s.mean('nScans')
        if fl is not None:
            fl.next('Raw E(p) - %s'%filename)
            fl.image(s)
        #{{{DC offset correction    
        s.ift('t2')
        if fl is not None:
            fl.next('raw E(p) time - %s'%filename)
            fl.image(s)
            #fl.show();quit()
        s *= exp(-abs((s.fromaxis('t2')-s.get_prop('acq_params')['tau_us']*1e-6))/Ep_apo)
        s.ft('t2')
        s.ift('t2')
        s.ift(['ph1'])
        Ep_t_max = s.getaxis('t2')[-1]
        rx_offset_corr = s['t2':(Ep_t_max*0.75,None)]
        rx_offset_corr = rx_offset_corr.mean(['t2'])
        s -= rx_offset_corr
        s.ft('t2')
        s.ft(['ph1'])
        zero_crossing = abs(select_pathway(s,Ep_signal_pathway)).C.sum('t2').argmin(powername,raw_index=True).item()
        #}}}
        if Ep_f_slice is None:
            for_lims = s
            drift_bounds = [
                    select_pathway(for_lims[powername,0],Ep_signal_pathway).argmax().item()
                    ,
                    select_pathway(s[powername,-1],Ep_signal_pathway).argmax().item()
                    ]
            drift_bounds.sort()
            this_Ep_f_slice = array(drift_bounds) + r_[-Ep_f_slice_auto_scale,Ep_f_slice_auto_scale]
        else:
            this_Ep_f_slice = Ep_f_slice
        s = s['t2':this_Ep_f_slice]
        s.ift('t2')
        #{{{phasing
        s.set_units('t2','s')
        best_shift,_ = hermitian_function_test(select_pathway(s,
            Ep_signal_pathway),
            echo_before=s.get_prop('acq_params')['tau_us']*1e-6*1.5,
            fl=fl)
        better = float("{:.6f}".format(best_shift))
        print(best_shift)
        actual_tau = s.get_prop('acq_params')['tau_us']/1e6
        if (best_shift < actual_tau-0.5e-3) or (best_shift > actual_tau+0.5e-3):
            if fl is not None:
                fl.text(r'\textcolor{red}{\textbf{I am hard-setting the first-order phase for dataset Ep}}')
                fl.basename =None
                best_shift = actual_tau
        print(best_shift)
        s.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
        s /= zeroth_order_ph(select_pathway(s['t2':0],Ep_signal_pathway))
        s.ft('t2')
        s.ift('t2')
        s = s['t2':(0,None)]
        s['t2':0] *= 0.5
        s.ft('t2')
        if fl is not None:
            fl.next('phased and FID sliced')
            fl.image(s)
            #fl.show();quit()
        #}}}
        #{{{Alignment
        last_p_max = select_pathway(s[powername,-1],Ep_signal_pathway).C.argmax('t2').item()
        first_p_max = select_pathway(s[powername,0],Ep_signal_pathway).C.argmax('t2').item()
        if (last_p_max <= 0) and (first_p_max <=0):
            last_p_max *= -1
            first_p_max *= -1
            if last_p_max > first_p_max:
                drift = last_p_max - first_p_max
            else:
                drift = first_p_max - last_p_max
        elif (last_p_max>=0) and (first_p_max >=0):
            if last_p_max > first_p_max:
                drift = last_p_max - first_p_max
            else:
                drift = first_p_max - last_p_max
        elif (last_p_max>=0) and(first_p_max <=0):
            first_p_max *= -1
            drift = first_p_max +last_p_max
        elif (last_p_max <=0) and(first_p_max >=0):
            last_p_max *= -1
            drift = last_p_max +first_p_max
        print("your Ep signal has a smear that spans %d Hz"%drift)    
        s[powername,:zero_crossing] *= -1
        if drift < Ep_drift_max:
            s.ift(['ph1'])
            opt_shift,sigma, my_mask = correl_align((s.C),
                    indirect_dim=powername,
                    signal_pathway=Ep_signal_pathway,
                    sigma = 1500)
            s.ift('t2')
            s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
            s.ft(['ph1'])
            s.ft('t2')
            s[powername,:zero_crossing] *= -1
        else:
            freq_diff = abs(select_pathway(s[powername,-1],
                Ep_signal_pathway)).argmax().item() - abs(select_pathway(
                    s[powername,0], Ep_signal_pathway)).argmax().item()
            freq_diff /= len(s.getaxis(powername))
            s.ift('t2')
            for i in range(len(s.getaxis(powername))):
                s[powername,i] *= np.exp(-1j*2*pi*(freq_diff+i*freq_diff)*s.fromaxis('t2'))
            s.ft('t2')
            ph0 = s['t2':(-900,0)].C.sum('t2')
            ph0 /= abs(ph0)
            s /= ph0
            s.ift(['ph1'])
            opt_shift,sigma, my_mask = correl_align((s),indirect_dim=powername,
                    signal_pathway=Ep_signal_pathway,sigma = 1500)
            s.ift('t2')
            s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
            s.ft('t2')
            s.ft(['ph1'])
        #}}}
        s.ift('t2')
        #{{{FID slice
        d=s.C
        d = d['t2':(0,None)]
        d['t2':0] *= 0.5
        d.ft('t2')
        if Ep_extra_zero:
            d.ift('t2')
            d /= zeroth_order_ph(select_pathway(d['t2':0],Ep_signal_pathway))
            d.ft('t2')
        if fl is not None:
            fl.next('E(p) FID slice after alignment')
            fl.image(d)
            #fl.show();quit()
        #}}}
        #{{{Integrate with error
        error_pathway = (set(((j) for j in range(ndshape(d)['ph1'])))
                - set(excluded_pathways)
                -set([(Ep_signal_pathway['ph1'])]))
        error_pathway = [{'ph1':j} for j in error_pathway]
        s_int,frq_slice = integral_w_errors(d.C,Ep_signal_pathway,error_pathway,
                cutoff=Ep_cutoff,
                indirect=powername,return_frq_slice=True)
        if fl is not None:
            fl.next('E(p) Integration Limits')
            fl.plot(select_pathway(d,Ep_signal_pathway))
            plt.axvline(x=frq_slice[0])
            plt.axvline(x=frq_slice[-1])
            #fl.show();quit()
        #}}}
        #{{{set time axis
        if log:
            s_int.getaxis(powername)[:]['start_times'] -= log_start_time
            s_int.getaxis(powername)[:]['stop_times'] -= log_start_time
            s_int.getaxis(powername)[-1]['stop_times'] =power_axis.getaxis(powername)[-1]
            time_axis = s_int.getaxis(powername)[:]['start_times']
            s_int.setaxis(powername,time_axis)
        else:
            s_int.setaxis(powername,[0.0, 0.28, 0.56, 0.89, 1.12, 1.41, 1.78, 2.0,
                2.24, 2.51, 2.82, 3.16, 3.55, 3.98])
        #}}}
        #{{{Normalize
        s_int /= s_int[powername,0:1]
        if Ep_flip:
            s_int[powername,1:] *= -1
        #}}}
        #{{{finding average power over steps
        if log:
            dnp_time_axis = s.C.getaxis(powername).copy()
            dnp_time_axis[0]['start_times'] = log_start_time
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
                if fl is not None:
                    fl.next('Instrument Power Log')
                    plt.axvline(x=time_start,color='red',alpha=0.5)
                    plt.axvline(x=time_stop,linestyle=':',color='red',alpha=0.5)
            avg_p_vs_t = nddata(power_vs_time.data,[-1],[powername])
            avg_p_vs_t.set_error(power_vs_time.get_error())
            avg_p_vs_t.setaxis(powername,dnp_time_axis[:]['start_times'])
            avg_p_vs_t.data[0] = 0
            avg_p_vs_t[powername,0].set_error(0)
            if fl is not None:
                fl.plot(avg_p_vs_t,'ro',capsize=6)
                plt.xlabel('Time / s')
                plt.ylabel('Power / W')
        #}}}
        #{{{set time axis to power for Ep
        s_int.rename(powername,'power')
        if log:
            power_axis = np.real(power_vs_time.data)
            s_int.setaxis('power',power_axis)
            s_int.set_error('power',power_vs_time.get_error())
        else:
            s_int.setaxis('power',hack_Ep_powers)
        enhancement=s_int
        if fl is not None:
            fl.next('Final Integrated Enhancement - %s'%filename)
            fl.plot(s_int['power',:-3],'o',capsize=6)
            fl.plot(s_int['power',-3:],'ro',capsize=6)
        #}}}
        #}}}
    if (has_IR and has_Ep):
        return T1p,enhancement
    elif (has_IR):
        return T1p
    elif (has_Ep):
        return enhancement
 