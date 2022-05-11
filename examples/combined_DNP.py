"""DNP Processing
=================

Take 2D dataset of enhancement and series of T1 datasets, and converts
to corresponding table of integrals. These are then combined to create
the $k_{\sigma}$s(p) curve for the given sample.
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
import sympy as sp
from sympy import symbols, Symbol
import os, time, h5py
import matplotlib.transforms as transforms
from Instruments.logobj import logobj
from itertools import cycle

init_logging(level="debug")

rcParams["image.aspect"] = "auto" #needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 5

color_cycle = cycle(['#1f77b4', '#ff7f0e', '#2ca02c',
    '#d62728', '#9467bd', '#8c564b', '#e377c2',
    '#7f7f7f', '#bcbd22', '#17becf'])
# {{{ this is just to show all the parameters
def list_symbs(f):
    list_symbs = []
    for j, k in f.output().items():
        s_repr = sp.latex(sp.Symbol(j))
        list_symbs.append(f"${s_repr} = {k:0.5g}$")
    list_symbs = "\n".join(list_symbs)
    return list_symbs
#}}}
fl = figlist_var()
#{{{variables to set
filename = "220429_50mM_TEMPOL_DNP_2"
file_location = "ODNP_NMR_comp/ODNP"
IR_sig_path = {'ph1':0,'ph2':1}
Ep_sig_path = {'ph1':1}
nPowers = 15
W = 0.14 + 1.024
R1w = 1/2.406166
dT1w = 0.2759
log_acquired = True
ppt = 1.516598e-3
#}}}
#{{{load in log if available
if log_acquired:
    with h5py.File(search_filename(filename+".h5", exp_type = file_location, unique = True),
            'r') as f:
        log_grp = f['log']
        thislog = logobj()
        thislog.__setstate__(log_grp)
        read_array = thislog.total_log
        read_dict = thislog.log_dict
        for j in range(len(read_array)):
            thistime, thisrx, thispower, thiscmd = read_array[j]
        fig, (ax_Rx,ax_power) = plt.subplots(2,1, figsize=(10,8))
        fl.next("Log Figure",fig=fig)
        ax_Rx.set_ylabel('Rx/mV')
        log_start_time = read_array['time'][0]
        relative_time = read_array['time']
        ax_Rx.plot(relative_time,read_array['power'],'.')
        mask = read_array['cmd'] != 0
        n_events = len(relative_time[mask])
        trans_power = transforms.blended_transform_factory(
                ax_power.transData, ax_power.transAxes)
        for j, thisevent in enumerate(read_array[mask]):
            ax_Rx.axvline(x=thisevent['time']-log_start_time)
            ax_power.axvline(x=thisevent['time']-log_start_time)
        ax_power.legend(**dict(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.))
        plt.tight_layout()
        power_axis = nddata(read_array['power'],[-1],['time'])
        power_axis.setaxis('time',relative_time-log_start_time)
        power_axis.name('power')
        fl.next('Instrument Power Log')
        coupler_atten = 22
        power_axis.data = (10**((power_axis.data+coupler_atten)/10+3))/1e6
        fl.plot(power_axis,'.')
#}}}
#{{{Process IR datasets
T1_list = []
powers = []
power_list = []
start_times = []
stop_times = []
errors = []
for (nodename, clock_correction, label) in [
        ('FIR_noPower', False, 'no_power'),
        ('FIR_30.0dBm', False, '30 dBm'),
        ('FIR_33.0dBm', False, '33 dBm'),
        ('FIR_34.5dBm', False, '34.5 dBm'),
        ('FIR_36.0dBm', False, '36.0 dBm')
        ]:
    fl.basename = "(%s)" % label
    thiscolor = next(color_cycle)
    IR = find_file(filename, 
            exp_type=file_location,
            expno = nodename,
            postproc = 'spincore_IR_v1',
            lookup=lookup_table)
    #{{{generate integrals
    IR_int, IR = generate_integrals(
            IR, signal_pathway = IR_sig_path,
            searchstr = label,
            f_range = (-0.25e3,0.15e3),
            indirect = 'vd',
            alias_slop=0,
            clock_correction = clock_correction,
            fl=None)
    #}}}
    #{{{Fit with lmfit
    M0, Mi, R1, vd = symbols("M_0 M_inf R_1 vd", real=True)
    functional_form = Mi*(1-(2-sp.exp(-W*R1))*sp.exp(-vd*R1))
    IR_data = nddata(IR_int.data,['vd'])
    IR_data.setaxis('vd',IR_int.getaxis('vd'))
    IR_fit = lmfitdata(IR_data)
    IR_fit.functional_form = functional_form
    IR_fit.set_guess(
            M_0 = dict(value=-5e6+1, min = -5e6, max = 0),
            M_inf = dict(value=5e6-1, min = 0, max = 5e6),
            R_1 = dict(value=30, min=0.001, max = 100)
            )
    fl.basename = None
    fl.next('IR fits - normalized')
    IR_fit.settoguess()
    guess = IR_fit.eval(100)
    IR_fit.fit()
    T1 = 1./IR_fit.output('R_1')
    Mi = IR_fit.output('M_inf')
    s_norm = IR_int/Mi
    fl.plot(s_norm, 'o', color=thiscolor,label = '%s normalized data'%label)
    T1_list.append(T1)
    fit = IR_fit.eval(100)
    fit /= Mi
    fl.plot(fit,
            color=thiscolor,
            alpha=0.5, label = '%s fit'%label)
    ax=plt.gca()
    plt.xlabel('vd/s')
    #}}}
    #{{{add power to power list
    if log_acquired:
        start_time = IR.get_prop('start_time')-log_start_time
        stop_time = IR.get_prop('stop_time')-log_start_time
        avg_power = power_axis['time':(start_time,stop_time)].mean('time',std=True)
        errors.append(avg_power.get_error())
        power_list.append(avg_power)
        powers.append(avg_power.data)
        fl.next('Instrument Power Log')
        start_times.append(start_time)
        stop_times.append(stop_time)
        plt.axvline(x=start_time,color='k')
        plt.axvline(x=stop_time,linestyle=':',color='k',alpha=0.5)
        nddata_p_vs_t = nddata(power_list,[-1], ['time'])
        nddata_p_vs_t.setaxis('time',start_times)
        for j in range(len(power_list)):
            nddata_p_vs_t['time',j] = power_list[j]
        nddata_p_vs_t.set_error(errors)
        fl.plot(nddata_p_vs_t, 'ko',capsize=6)
    else:
        power_dB = IR.get_prop('acq_params')['meter_power']
        power_W = (10**((power_dB+coupler_atten)/10+3))/1e6
        power_list.append(power_W)
    #}}}
fl.show();quit()    
#{{{ process enhancement dataset
for nodename,postproc in [
        ('enhancement','spincore_ODNP_v3')
        ]:
    Ep = find_file(filename, exp_type=file_location, expno=nodename,
            postproc = postproc,lookup=lookup_table)
    Ep.reorder(['ph1','time'])
    #{{{set power axis based on log
    if log_acquired:
        dnp_t_axis = Ep.C.getaxis('time').copy()
        dnp_t_axis[0]['start_times'] = log_start_time
        dnp_t_axis[:]['start_times'] -= log_start_time
        dnp_t_axis[:]['stop_times'] -= log_start_time
        nddata_time_axis = nddata(dnp_t_axis,[-1],['time'])
        new_t_axis = nddata_time_axis.C.data
        new_t_axis = nddata(new_t_axis,[-1],['time'])
        power_vs_time = ndshape(nddata_time_axis).alloc().set_units('time','s')
        power_vs_time.set_error(0)
        power_vs_time.setaxis('time',new_t_axis.data)
        fl.next('Instrument Power Log')
        for j, (time_start,time_stop) in enumerate(zip(dnp_t_axis[:]['start_times'],dnp_t_axis[:]['stop_times'])):
            power_vs_time['time',j] = power_axis['time':((time_start),(time_stop))].mean('time',std=True)
            power_vs_time.set_units('time','s')
            plt.axvline(x=time_start,color='red',alpha=0.5)
            plt.axvline(x = time_stop, linestyle=':', color='red',alpha=0.5)
        avg_p_vs_t = nddata(power_vs_time.data, [-1],['time'])
        avg_p_vs_t.set_error(power_vs_time.get_error())
        avg_p_vs_t.setaxis('time',dnp_t_axis[:]['start_times'])
        avg_p_vs_t.data[0] = 0
        avg_p_vs_t['time',0].set_error(0)
        fl.plot(avg_p_vs_t,'ro',capsize=6)
        plt.xlabel('Time / s')
        plt.ylabel('Power / W')
        Ep.rename('time','power')
        power_axis = np.real(power_vs_time.data)
        Ep.setaxis('power',power_axis)
        Ep.set_error('power',power_vs_time.get_error())
    #}}}
    #{{{set power axis based on meter readings
    else:
        power_axis_dBm = array(Ep.get_prop('acq_params')['power_settings_dBm'])
        power_axis_W = zeros_like(power_axis_dBm)
        power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
        power_axis_W = r_[0,power_axis_W]
        Ep.setaxis('power',power_axis_W)
    #}}}    
    Ep_int,Ep = generate_integrals(Ep,signal_pathway = Ep_sig_path,
            searchstr = 'E(p)',
            f_range = (-500,500),
            indirect = 'power',
            alias_slop = 0,
            clock_correction = False,
            fl=fl)
#}}}
#{{{make k_rho_inv and R1p
T10 = T1_list[0]
R10 = T10**-1
T1_list = np.array(T1_list)
T1s = nddata(T1_list,['power'])
T1s.setaxis('power', np.array(powers))
R1p = T1s ** -1
p = np.array(powers)
T1w = R1w **-1
fl.next(r'$T_{1}$(p)')
fl.plot(T1s,'o')
polyorder = 4
R1water = 1 / (T1w + p * dT1w)
rho_inv = (1.0/T1_list-R1water)**-1
rho_inv_nddata = nddata(rho_inv,[-1],['power'])
rho_inv_nddata.setaxis('power',p)
#{{{shift axis so polyorder is centered
shift_value = rho_inv_nddata.C.getaxis('power')[2]
rho_inv_nddata.setaxis('power',lambda x: x-shift_value)
pfine = nddata(np.linspace(rho_inv_nddata.getaxis('power')[0],rho_inv_nddata.getaxis('power')[-1],nPowers),'power')
coeff = rho_inv_nddata.polyfit('power',order=polyorder)
p_fine = nddata(np.linspace(0,p[-1],nPowers),'power')
T1water_fine = np.array(T1w + p_fine*dT1w)
R1water_fine = T1water_fine ** -1
rho_inv_fine = 0
for j in range(polyorder + 1):
    rho_inv_fine += coeff[j] * pfine**j
rho_inv_fine.setaxis('power',lambda x: x+shift_value)
rho_inv_nddata.setaxis('power',lambda x: x+shift_value)
fl.next(r'$(k_{\rho})^{-1}$')
fl.plot(rho_inv_nddata,'o')
fl.plot(rho_inv_fine)
R1p_fine = rho_inv_fine**-1 + R1water_fine
fl.next(r'$R_{1}$(p)')
fl.plot(R1p,'o')
fl.plot(R1p_fine)
ksigs_T = (ppt)*(1-Ep_int['power',:-3])*(R1p_fine)
ksC, phalf, power = symbols("ksC phalf power", real =True)
functional_form = ksC*power/(phalf+power)
fl.next(r'$k_{\sigma}*s(p)*C$')
fl.plot(ksigs_T,'o', alpha=0.5, label='data')
this_data = nddata(ksigs_T.real.data,['power'])
this_data.setaxis('power',ksigs_T.getaxis('power'))
sigma_fit = lmfitdata(this_data)
sigma_fit.functional_form=functional_form
sigma_fit.set_guess(
        phalf=dict(value=0.5, mn=0, max = 5),
        ksC=dict(value=1,min=0,max =500)
        )
sigma_fit.settoguess()
guess_line = sigma_fit.eval(100)
sigma_fit.fit()
fit_line = fl.plot(sigma_fit.eval(100),
        alpha=0.5, label = 'lmfit fit')
ax = plt.gca()
text(0.5,0.5, "Result: %s" % sigma_fit.latex(),
        ha = 'center',
        va = 'center',
        transform = ax.transAxes
        )
text(0.5,0.5,(3 *"\n") + list_symbs(sigma_fit),
        ha = 'center',
        va = 'top',
        size = 12,
        transform = ax.transAxes
        )
fl.show();quit()
