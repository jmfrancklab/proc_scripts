import os, time, h5py
import pylab as plt
from numpy import empty
from sympy import exp as s_exp
from sympy import symbols, Symbol, latex
from matplotlib.ticker import FuncFormatter
import matplotlib.transforms as transforms
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from Instruments.logobj import logobj
from scipy.io import loadmat
from itertools import cycle
fl = figlist_var()
color_cycle = cycle(['red','orange','yellow','green','cyan','blue','purple','magenta',
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
    '#7f7f7f', '#bcbd22', '#17becf'])

def generate_T1_Ep(filename,
        has_Ep = True,
        has_IR = True,
        Ep_signal_pathway={'ph1':1},
        IR_signal_pathway = {'ph1':0,'ph2':1},
        Ep_nodename = 'enhancement',
        IR_nodenames = [('FIR_noPower'),('FIR_30.0dBm'),('FIR_33.0dBm'),('FIR_34.5dBm'),('FIR_36.0dBm')],
        Ep_f_slice = (-500,500),
        IR_f_slice = (-500,500),
        Mi_max = 1e6,
        drift_max = 100,
        excluded_pathways = [(0,0)],
        Ep_postproc = 'spincore_ODNP_v3',
        IR_postproc='spincore_IR_v1',
        nPowers = 15,
        W = 7+1.024,
        IR_cutoff = 0.15,
        Ep_cutoff = 0.15,
        log=True,
        Ep_flip = False, 
        clock_correction = True,
        older = False,
        IR_full_flip=False,
        push_zero_time = 0,
        coupler_atten = 22,
        fl=None):
    T1_list = []
    errors = []
    power_list = []
    start_times = []
    powers = []
    if log:
    #{{{load in log
        my_filename = search_filename(filename+".h5",exp_type='ODNP_NMR_comp/ODNP',unique=True)
        with h5py.File(my_filename,'r') as f:
            #logging.debug(strm("for log, loading file",my_filename))
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
        
        #trans_power = transforms.blended_transform_factory(
        #        ax_power.transData,ax_power.transAxes)
        #trans_Rx = transforms.blended_transform_factory(
        #        ax_Rx.transData, ax_Rx.transAxes)
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
    if has_IR:
    #{{{IR processing
        for (j,nodename) in enumerate(IR_nodenames):
            IR = find_file(filename, exp_type='ODNP_NMR_comp/ODNP', expno=nodename,
                    postproc=IR_postproc, lookup=lookup_table)
            if older:
                IR = IR['indirect',j]
            else:
                pass
            IR.reorder(['ph1','ph2','vd'])
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
            if IR_f_slice == None:
                for_lims = IR.C.mean('nScans')
                right_f = select_pathway(for_lims['vd',0],IR_signal_pathway).argmax() + 200
                left_f = select_pathway(IR.C.mean('nScans')['vd',-1],IR_signal_pathway).argmax() - 200
                IR_f_slice = (left_f.data,right_f.data)
            zero_crossing = abs(select_pathway(IR['t2':IR_f_slice].C.mean('nScans'),IR_signal_pathway)).C.sum('t2').argmin('vd',raw_index=True).item()
            if fl is not None:
                fl.next('Raw IR for %s'%nodename)
                fl.image(IR.C.mean('nScans'))
            IR=IR['t2':IR_f_slice]
            IR.ift('t2')
            #{{{clock correction
            clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
            IR.ft('t2')
            if (len(IR.getaxis('nScans'))) > 1:
                s_clock = select_pathway(IR.C.mean('nScans'),IR_signal_pathway).C.sum('t2')
            else:
                s_clock = select_pathway(IR.C, IR_signal_pathway).C.sum('t2')
            IR.ift(['ph1','ph2'])
            min_index = abs(s_clock).argmin('vd',raw_index=True).item()
            s_clock *= np.exp(-1j*clock_corr*IR.fromaxis('vd'))
            s_clock['vd',:min_index+1] *=-1
            s_clock.sum('vd').run(abs)
            clock_corr = s_clock.argmax('clock_corr').item()
            IR *= np.exp(-1j*clock_corr*IR.fromaxis('vd'))
            IR.ft(['ph1','ph2'])
            IR.ift('t2')
            #}}}
            #{{{phasing
            best_shift = hermitian_function_test(select_pathway(IR,
                IR_signal_pathway),aliasing_slop=0)
            actual_tau = IR.get_prop('acq_params')['tau_us']/1e6
            if (best_shift < actual_tau-1e-3) or (best_shift > actual_tau+1e-3):
                best_shift = actual_tau
                if fl is not None:
                    fl.text(r'\textcolot{red}{\textbf{I am hard-setting the first-order phase for dataset %s}}'%nodename)
            IR.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
            IR /= zeroth_order_ph(select_pathway(IR['t2':0],IR_signal_pathway))
            IR.ft('t2')
            #}}}
            #{{{Alignment
            last_vd_max = select_pathway(IR['vd',-1].C.mean('nScans'),IR_signal_pathway).C.argmax('t2').item()
            first_vd_max = select_pathway(IR['vd',0].C.mean('nScans'),IR_signal_pathway).C.argmax('t2').item()
            if (last_vd_max <= 0) and (first_vd_max <= 0):
                last_vd_max *= -1
                first_vd_max *= -1
                if last_vd_max > first_vd_max:
                    drift = last_vd_max - first_vd_max
                else:
                    drift = first_vd_max - last_vd_max
            elif (last_vd_max>=0) and (first_vd_max >=0):
                if last_vd_max > first_vd_max:
                    drift = last_vd_max - first_vd_max
                else:
                    drift = first_vd_max - last_vd_max
            elif (last_vd_max>=0) and(first_vd_max <=0):
                first_vd_max *= -1
                drift = first_vd_max +last_vd_max
            elif (last_vd_max <=0) and(first_vd_max >=0):
                last_vd_max *= -1
                drift = last_vd_max +first_vd_max
            print("your signal has a smear that spans %d Hz"%drift)
            mysgn = determine_sign(select_pathway(IR.C,IR_signal_pathway))
            if drift < drift_max:
                IR.ift(['ph1','ph2'])
                opt_shift,sigma, my_mask = correl_align((IR.C*mysgn),indirect_dim='vd',
                        signal_pathway=IR_signal_pathway,sigma=1200)#again, bad data
                IR.ift('t2')
                IR *= np.exp(-1j*2*pi*opt_shift*IR.fromaxis('t2'))
                IR.ft(['ph1','ph2'])
                IR.ft('t2')
            else:
                IR['vd',:zero_crossing] *= -1
                freq_diff = abs(select_pathway(IR['nScans',0]['vd',-1],IR_signal_pathway)).argmax().item() - abs(select_pathway(IR['nScans',0]['vd',0], IR_signal_pathway)).argmax().item()
                freq_diff /= len(IR.getaxis('vd'))
                IR.ift('t2')
                for i in range(len(IR.getaxis('vd'))):
                    IR['vd',i] *= np.exp(-1j*2*pi*(freq_diff+i*freq_diff)*IR.fromaxis('t2'))
                IR.ft('t2')
                IR['vd',:zero_crossing] *= -1
                ph0 = IR['t2':(-900,0)].C.sum('t2')
                ph0 /= abs(ph0)
                IR /= ph0
                IR.ift(['ph1','ph2'])
                opt_shift,sigma, my_mask = correl_align((IR),indirect_dim='vd',
                        signal_pathway=IR_signal_pathway,sigma = 500)
                IR.ift('t2')
                IR *= np.exp(-1j*2*pi*opt_shift*IR.fromaxis('t2'))
                IR.ft('t2')
                IR.ft(['ph1','ph2'])
                IR['vd',zero_crossing:] *= -1
            #}}}
            IR.ift('t2')
            #{{{FID slice
            d=IR.C
            d = d['t2':(0,None)]
            d['t2':0] *= 0.5
            d.ft('t2')
            if fl is not None:
                fl.next('IR aligned and FID sliced - %s'%nodename)
                fl.image(d)
            #}}}
            #{{{Integrate with error
            this_IR = (d.C*mysgn).mean('nScans')
            error_pathway = (set(((j,k) for j in range(ndshape(d)['ph1']) for k in range(ndshape(d)['ph2'])))
                    - set(excluded_pathways)
                    -set([(IR_signal_pathway['ph1'],IR_signal_pathway['ph2'])]))
            error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
            s_int,frq_slice = integral_w_errors(this_IR,IR_signal_pathway,error_pathway,
                    cutoff = IR_cutoff,
                    indirect='vd',return_frq_slice=True)
            if fl is not None:
                fl.next('IR Integration Limits - %s'%nodename)
                for j in range(len(d.getaxis('vd'))):
                    fl.plot(select_pathway((d.C*mysgn).mean('nScans')['vd',j],
                        IR_signal_pathway),label = '%d'%j)
                plt.axvline(frq_slice[0])
                plt.axvline(frq_slice[-1])
            #}}}
            #{{{Fitting
            s_int['vd',:zero_crossing] *= -1
            thiscolor = next(color_cycle)
            if IR_full_flip:
                if s_int['vd',-1] <0:
                    s_int *= -1
            M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd",real=True)
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
            if fl is not None:
                fl.next('IR fit - before norm')
                fl.plot(s_int, 'o', color=thiscolor, label='%s'%nodename)
                fl.plot(fit, ls = this_ls, color=thiscolor,alpha=0.5)
                fl.next('IR fit - normalized')
                fl.plot(s_int/Mi, 'o', color=thiscolor, label=nodename)
                fl.plot(fit/Mi, 
                        ls=this_ls,
                        color=thiscolor,
                        alpha=0.5,
                        label='fit for %s'%nodename,
                        )
                ax=plt.gca()
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
            else:
                power_dB = IR.get_prop('acq_params')['meter_power']
                power_W = (1e-2*10**((power_dB)*1e-1))*1e-1 
                powers.append(power_W)
                power_list.append(power_W)
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
    #}}}
    #}}}
    if has_Ep:
    #{{{Ep processing    
        s = find_file(filename, exp_type = 'ODNP_NMR_comp/ODNP', expno= Ep_nodename, 
                postproc=Ep_postproc, lookup=lookup_table)
        s.reorder(['ph1','time'])
        if fl is not None:
            fl.next('Raw E(p)')
            fl.image(s)
        #{{{DC offset correction    
        s.ift('t2')
        s.ift(['ph1'])
        Ep_t_max = s.getaxis('t2')[-1]
        rx_offset_corr = s['t2':(Ep_t_max*0.75,None)]
        rx_offset_corr = rx_offset_corr.mean(['t2'])
        s -= rx_offset_corr
        s.ft('t2')
        s.ft(['ph1'])
        #}}}
        s = s['t2':Ep_f_slice]
        s.ift('t2')
        #{{{phasing
        best_shift = hermitian_function_test(select_pathway(s,
            Ep_signal_pathway),aliasing_slop=0)
        actual_tau = s.get_prop('acq_params')['tau_us']/1e6
        if (best_shift < actual_tau-1e-3) or (best_shift > actual_tau+1e-3):
            best_shift = actual_tau
            if fl is not None:
                fl.text(r'\textcolor{red}{\textbf{I am hard-setting the first-order phase for dataset %s}}'%datasetname)
        s.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
        s /= zeroth_order_ph(select_pathway(s['t2':0],Ep_signal_pathway))
        s.ft('t2')
        #}}}
        #{{{Alignment
        last_p_max = select_pathway(s['time',-1],Ep_signal_pathway).C.argmax('t2').item()
        first_p_max = select_pathway(s['time',0],Ep_signal_pathway).C.argmax('t2').item()
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
        mysgn = determine_sign(select_pathway(s.C,Ep_signal_pathway))
        if drift < drift_max:
            s.ift(['ph1'])
            opt_shift,sigma, my_mask = correl_align((s.C*mysgn),
                    indirect_dim='time',
                    signal_pathway=Ep_signal_pathway)
            s.ift('t2')
            s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
            s.ft(['ph1'])
            s.ft('t2')
        else:
            freq_diff = abs(select_pathway(s['nScans',0]['time',-1],
                Ep_signal_pathway)).argmax().item() - abs(select_pathway(
                    s['nScans',0]['time',0], Ep_signal_pathway)).argmax().item()
            freq_diff /= len(s.getaxis('time'))
            s.ift('t2')
            for i in range(len(s.getaxis('time'))):
                s['time',i] *= np.exp(-1j*2*pi*(freq_diff+i*freq_diff)*s.fromaxis('t2'))
            s.ft('t2')
            ph0 = s['t2':(-900,0)].C.sum('t2')
            ph0 /= abs(ph0)
            s /= ph0
            s.ift(['ph1'])
            opt_shift,sigma, my_mask = correl_align((s),indirect_dim='time',
                    signal_pathway=Ep_signal_pathway,sigma = 150)
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
        if fl is not None:
            fl.next('FID slice after alignment')
            fl.image(d)
        #}}}
        #{{{Integrate with error
        error_pathway = (set(((j) for j in range(ndshape(d)['ph1'])))
                - set(excluded_pathways)
                -set([(Ep_signal_pathway['ph1'])]))
        error_pathway = [{'ph1':j} for j in error_pathway]
        s_int,frq_slice = integral_w_errors(d.C.mean('nScans')*mysgn.mean('nScans'),Ep_signal_pathway,error_pathway,
                cutoff=Ep_cutoff,
                indirect='time',return_frq_slice=True)
        if fl is not None:
            fl.next('E(p) Integration Limits')
            fl.plot(select_pathway(d.C.mean('nScans'),Ep_signal_pathway))
            plt.axvline(x=frq_slice[0])
            plt.axvline(x=frq_slice[-1])
        #}}}
        #{{{set time axis
        if log:
            s_int.getaxis('time')[:]['start_times'] -= log_start_time
            s_int.getaxis('time')[:]['stop_times'] -= log_start_time
            s_int.getaxis('time')[-1]['stop_times'] =power_axis.getaxis('time')[-1]
            time_axis = s_int.getaxis('time')[:]['start_times']
            s_int.setaxis('time',time_axis)
        else:
            if ('power_settings') in s.get_prop('acq_params'):
                power_dB = s.get_prop('acq_params')['power_settings']
                powers_W = [0]
            else:
                power_dB = s.get_prop('acq_params')['power_settings_dBm']
            for j in range(len(power_dB)):
                power_W = (1e-2*10**((power_dB[j])+10.)*1e-1)
                powers_W.append(power_W)
            s_int.setaxis('time',power_W)
        #}}}
        #{{{Normalize
        s_int['time',:] /= s_int.data[0]
        if Ep_flip:
            s_int['time',1:] *= -1
        #}}}
        #{{{finding average power over steps
        if log:
            dnp_time_axis = s.C.getaxis('time').copy()
            dnp_time_axis[0]['start_times'] = log_start_time
            dnp_time_axis[:]['start_times'] -= log_start_time
            dnp_time_axis[:]['stop_times'] -= log_start_time
            nddata_time_axis = nddata(dnp_time_axis,[-1],['time'])
            new_time_axis = nddata_time_axis.C.data
            new_time_axis = nddata(new_time_axis,[-1],['time'])
            power_vs_time = ndshape(nddata_time_axis).alloc().set_units('time','s')
            power_vs_time.set_error(0)
            power_vs_time.setaxis('time',new_time_axis.data)
        #}}}
        #{{{find values for Ep
            for j,(time_start,time_stop) in enumerate(zip(dnp_time_axis[:]['start_times'],dnp_time_axis[:]['stop_times'])):
                power_vs_time['time',j] = power_axis['time':((time_start),(time_stop))].mean('time',std=True)
                power_vs_time.set_units('time','s')
                if fl is not None:
                    fl.next('Instrument Power Log')
                    plt.axvline(x=time_start,color='red',alpha=0.5)
                    plt.axvline(x=time_stop,linestyle=':',color='red',alpha=0.5)
            avg_p_vs_t = nddata(power_vs_time.data,[-1],['time'])
            avg_p_vs_t.set_error(power_vs_time.get_error())
            avg_p_vs_t.setaxis('time',dnp_time_axis[:]['start_times'])
            avg_p_vs_t.data[0] = 0
            avg_p_vs_t['time',0].set_error(0)
            if fl is not None:
                fl.plot(avg_p_vs_t,'ro',capsize=6)
                plt.xlabel('Time / s')
                plt.ylabel('Power / W')
        #}}}
        #{{{set time axis to power for Ep
        s_int.rename('time','power')
        if log:
            power_axis = np.real(power_vs_time.data)
            s_int.setaxis('power',power_axis)
            s_int.set_error('power',power_vs_time.get_error())
        else:
            pass
        enhancement=s_int
        if fl is not None:
            fl.next('Final Integrated Enhancement')
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
 
