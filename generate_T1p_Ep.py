import pylab as plt
from sympy import exp as s_exp
from sympy import symbols
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from itertools import cycle
fl = figlist_var()
color_cycle = cycle(['#1f77b4', '#ff7f0e', '#2ca02c',
    '#d62728', '#9467bd', '#8c564b', '#e377c2',
    '#7f7f7f', '#bcbd22', '#17becf'])

logger = init_logging("info")
Ep_signal_pathway = {'ph1':1}
IR_signal_pathway = {'ph1':0,'ph2':-1}
fl=fl_mod()
filename='220616_ras_I21_capProbe_ODNP'
file_location = 'ODNP_NMR_comp/ODNP'
Ep_f_slice=(17,317)
excluded_pathways = [(0,0)]
W = 4+1.024
h5filename = "processed_ras.h5"
conctag = "I21_220616"
T1_list = []
power_list = []
start_times = []
stop_times = []
errors=[]
FIR=True
Ep = True
#{{{make T1(p)
if FIR:
    for nodename, postproc, IR_f_slice, power in [
            ('FIR_noPower','spincore_IR_v1',(-47,253),0),
            ('FIR_32dBm','spincore_IR_v1',(-11,289),1.58),
            ('FIR_36dBm','spincore_IR_v1',(-5,130),3.98),
            ]:
        IR = find_file(filename,exp_type=file_location,expno=nodename,
                postproc=postproc,lookup=lookup_table)
        IR.reorder(['ph1','ph2','vd'])
        fl.next('Raw IR')
        fl.image(IR)
        #fl.show();quit()#########################################################
        IR['ph2',0]['ph1',0]['t2':0]=0 #kill axial noise
        IR.ift('t2')
        IR.ift(['ph1','ph2'])
        t_max = IR.C.getaxis('t2')[-1]
        rx_offset_corr = IR['t2':(t_max*0.75,None)]
        rx_offset_corr = rx_offset_corr.mean(['t2'])
        IR -= rx_offset_corr
        IR.ft('t2')
        IR.ft(['ph1','ph2'])
        zero_crossing = abs(select_pathway(IR['t2':IR_f_slice].C.mean('nScans'),IR_signal_pathway)).C.sum('t2').argmin('vd',raw_index=True).item()
        IR=IR['t2':IR_f_slice]
        IR.ift('t2')
        #{{{clock correction
        clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
        IR.ft('t2')
        fl.next('before clock correction')
        fl.image(IR.C.setaxis('vd','#').set_units('vd','scan #'))
        s_clock=select_pathway(IR.C,IR_signal_pathway).sum('t2')
        IR.ift(['ph1','ph2'])
        min_index = abs(s_clock).argmin('vd',raw_index=True).item()
        s_clock *= np.exp(-1j*clock_corr*IR.fromaxis('vd'))
        s_clock['vd',:min_index+1] *=-1
        s_clock.sum('vd').run(abs)
        fl.next('clock correction')
        fl.plot(s_clock,'.',alpha=0.7)
        clock_corr = s_clock.argmax('clock_corr').item()
        plt.axvline(x=clock_corr, alpha=0.5, color='r')
        IR *= np.exp(-1j*clock_corr*IR.fromaxis('vd'))
        IR.ft(['ph1','ph2'])
        fl.next('after auto-clock correction')
        fl.image(IR.C.setaxis('vd','#'))
        IR.ift('t2')
        #}}}
        #{{{phasing
        IR.ft('t2')
        IR.ift('t2')
        best_shift = hermitian_function_test(select_pathway(IR,IR_signal_pathway),aliasing_slop=0)
        print("estimated shift for IR is", best_shift)
        actual_tau = IR.get_prop('acq_params')['tau_us']/1e6
        if (best_shift < actual_tau-1e-3) or (best_shift > actual_tau+1e-3):
            print("go back and try a different slice. For now I am setting best shift = actual tau")
            best_shift = actual_tau
        IR.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
        IR.ft('t2')
        IR /= zeroth_order_ph(select_pathway(IR['t2':0],IR_signal_pathway))
        fl.next('phased')
        fl.image(IR)
        #fl.show();quit()####################################################################
        IR.ift('t2')
        #}}}
        #{{{Alignment
        IR.ft('t2')
        last_vd_max = select_pathway(IR['vd',-1],IR_signal_pathway).C.argmax('t2').item()
        first_vd_max = select_pathway(IR['vd',0],IR_signal_pathway).C.argmax('t2').item()
        if (last_vd_max < 0) and (first_vd_max <0):
            last_vd_max *= -1
            first_vd_max *= -1
            if last_vd_max > first_vd_max:
                drift = last_vd_max - first_vd_max
            else:
                drift = first_vd_max - last_vd_max
        elif (last_vd_max>0) and (first_vd_max >0):
            if last_vd_max > first_vd_max:
                drift = last_vd_max - first_vd_max
            else:
                drift = first_vd_max - last_vd_max
        elif (last_vd_max>0) and(first_vd_max <0):
            first_vd_max *= -1
            drift = first_vd_max +last_vd_max
        elif (last_vd_max <0) and(first_vd_max >0):
            last_vd_max *= -1
            drift = last_vd_max +first_vd_max
print("your signal has a smear that spans %d Hz"%drift)    
        if drift < 60:
            mysgn = determine_sign(select_pathway(IR,IR_signal_pathway))
            IR.ift(['ph1','ph2'])
            opt_shift,sigma, my_mask = correl_align((IR*mysgn),indirect_dim='vd',
                    signal_pathway=IR_signal_pathway,sigma=1500)
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
            IR['vd',:zero_crossing] *= -1
        fl.next('Aligned')
        fl.image(IR)
        #fl.show();quit()################################################################ 
        #}}}
        IR.ift('t2')
        #{{{FID slice
        d=IR.C
        d = d['t2':(0,None)]
        d['t2':0] *= 0.5
        d.ft('t2')
        fl.next('FID slice')
        fl.image(d)
        d.ift('t2')
        d /= zeroth_order_ph(select_pathway(d['t2':0],IR_signal_pathway))
        d.ft('t2')
        #fl.show();quit()
        #}}}
        #{{{Integrate with error
        error_pathway = (set(((j,k) for j in range(ndshape(d)['ph1']) for k in range(ndshape(d)['ph2'])))
                - set(excluded_pathways)
                -set([(IR_signal_pathway['ph1'],IR_signal_pathway['ph2'])]))
        error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
        s_int,frq_slice = integral_w_errors(d.C.mean('nScans'),IR_signal_pathway,error_pathway,
                cutoff = 0.1,
                indirect='vd',return_frq_slice=True)
        fl.next('diag')
        fl.plot(select_pathway(d.C.mean('nScans'),IR_signal_pathway))
        plt.axvline(frq_slice[0])
        plt.axvline(frq_slice[-1])
        #fl.show();quit()
        #}}}
        #{{{Fitting Routine
        M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd",real=True)
        functional_form = Mi*(1-(2-s_exp(-W*R1))*s_exp(-vd*R1))
        IR_data = nddata(s_int.data,['vd'])
        IR_data.setaxis('vd',s_int.getaxis('vd'))
        f = lmfitdata(IR_data)
        f.functional_form = functional_form
        f.set_guess(
                M_0 = dict(value=-1e6+1, min = -1e6, max = 0),
                M_inf = dict(value=1e6-1, min = 0, max=1e6),
                R_1 = dict(value=1.0, min=0.001, max=100)
                )
        fl.next('IR fit')
        fl.plot(s_int, 'o', label='data for %s'%nodename)
        f.settoguess()
        guess = f.eval(100)
        f.fit()
        T1 = 1./f.output('R_1')
        T1_list.append(T1)
        power_list.append(power)
        this_ls = "-"
        fit_line = fl.plot(f.eval(100), 
                ls=this_ls,
                alpha=0.5,
                label='fit for %s'%nodename,
                )
        ax=plt.gca()
        print("T1 for %s: %0.4f"%(nodename, T1))
        #fl.show();quit()
        #}}}
    #{{{finding average power for each T1
    logger.info(strm("T1 list:",T1_list)) 
    logger.info(strm("Power list:",power_list))
    T10 = T1_list[0]
    T1p = nddata(T1_list,[-1],['power'])
    T1p.setaxis('power',power_list)
    fl.next(r'$T_{1}$(p) vs power')
    fl.plot(T1p.data,'ro')
    fl.show();quit()##################################################################3333
    #T1p.name('T1p').hdf5_write(f'{h5filename}/{conctag}')
    #}}}
#}}}    
#{{{Load process enhancement
if Ep:
    for nodename,postproc in [
            ('enhancement','spincore_ODNP_v1')
            ]:
        #{{{
        s = find_file(filename,exp_type=file_location, expno=nodename,
                postproc=postproc,lookup=lookup_table)
        #{{{for old processing
        s.reorder(['ph1','power'])
        fl.next('Raw Ep')
        fl.image(s)
        #fl.show();quit()############################################################
        s.ift('t2')
        s.ift(['ph1'])
        Ep_t_max = s.C.getaxis('t2')[-1]
        rx_offset_corr = s['t2':(Ep_t_max*0.75,None)]
        rx_offset_corr = rx_offset_corr.mean(['t2'])
        s -= rx_offset_corr
        s.ft('t2')
        s.ft(['ph1'])
        s = s['t2':Ep_f_slice]
        s.ift('t2')
        #}}}
        #{{{phasing
        best_shift = hermitian_function_test(select_pathway(s,
            Ep_signal_pathway),aliasing_slop=0)
        print("estimated best_shift:",best_shift)
        actual_tau = s.get_prop('acq_params')['tau_us']/1e6
        if (best_shift < actual_tau-1e-3) or (best_shift > actual_tau+1e-3):
            print("go back and try a different slice. For now I am setting best shift = actual tau")
            best_shift = actual_tau
        s.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
        s.ft('t2')
        s /= zeroth_order_ph(select_pathway(s,Ep_signal_pathway))
        s.ift('t2')
        s.ft('t2')
        fl.next('Ep phase corrected')
        fl.image(s)
        #fl.show();quit()###################################################################
        #}}}
        #{{{Alignment
        last_p_max = select_pathway(s['power',-1],Ep_signal_pathway).C.argmax('t2').item()
        first_p_max = select_pathway(s['power',0],Ep_signal_pathway).C.argmax('t2').item()
        if first_p_max < 0:
            first_span = 0-(first_p_max* -1)
            drift = last_p_max + first_span
        else:
            if last_p_max > first_p_max:
                drift = last_p_max - first_p_max
            else:
                drift = first_p_max - last_p_max
        print("your signal has a smear that spans %d Hz"%drift)    
        if drift < 200:
            mysgn = determine_sign(select_pathway(s,Ep_signal_pathway))
            s.ift(['ph1'])
            opt_shift,sigma, my_mask = correl_align((s*mysgn),indirect_dim='power',
                    signal_pathway=Ep_signal_pathway,sigma=1500)
            s.ift('t2')
            s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
            s.ft(['ph1'])
            s.ft('t2')
        else:
            freq_diff = abs(select_pathway(s['nScans',0]['power',-1],
                Ep_signal_pathway)).argmax().item() - abs(select_pathway(
                    s['nScans',0]['power',0], Ep_signal_pathway)).argmax().item()
            freq_diff /= len(s.getaxis('power'))
            s.ift('t2')
            for i in range(len(s.getaxis('power'))):
                s['power',i] *= np.exp(-1j*2*pi*(freq_diff+i*freq_diff)*s.fromaxis('t2'))
            s.ft('t2')
            ph0 = s['t2':(-900,0)].C.sum('t2')
            ph0 /= abs(ph0)
            s /= ph0
            s.ift(['ph1'])
            opt_shift,sigma, my_mask = correl_align((s),indirect_dim='power',
                    signal_pathway=Ep_signal_pathway,sigma = 1500)
            s.ift('t2')
            s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
            s.ft('t2')
            s.ft(['ph1'])
        fl.next('Aligned')
        fl.image(s)
        #fl.show();quit()################################################################ 
        #}}}
        s.ift('t2')
        #{{{FID slice
        d=s.C
        d = d['t2':(0,None)]
        d['t2':0] *= 0.5
        d.ft('t2')
        #}}}
        #{{{Integrate with error
        error_pathway = (set(((j) for j in range(ndshape(d)['ph1'])))
                - set(excluded_pathways)
                -set([(Ep_signal_pathway['ph1'])]))
        error_pathway = [{'ph1':j} for j in error_pathway]
        s_int,frq_slice = integral_w_errors(d.C.mean('nScans'),Ep_signal_pathway,error_pathway,
                indirect='power',return_frq_slice=True)
        fl.next('Ep diag')
        fl.plot(select_pathway(d.C.mean('nScans'),Ep_signal_pathway))
        plt.axvline(x=frq_slice[0])
        plt.axvline(x=frq_slice[-1])
        #fl.show();quit()################################################################
        #}}}
        #{{{Normalize and flip
        s_int['power',:] /= s_int.data[0]
        fl.next('E(p) before power correction')
        fl.plot(s_int['power',:-3],'ko',capsize=2,alpha=0.3)
        fl.plot(s_int['power',-3:],'ro',capsize=2,alpha=0.3)
        fl.show();quit()###############################################################3
        #}}}
        #}}}
        s_int.name('Ep').hdf5_write(f'{h5filename}/{conctag}')
        #}}}


