from pylab import *
from pyspecdata import *
from pyspecProcScripts import hermitian_function_test, zeroth_order_ph
from pyspecProcScripts import integrate_limits, integral_w_errors
fl = figlist_var()

for filename, expno, subplot_name, in [
        #('210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_1000_before_noApod_ex','prealignment','tau1000, before'),
        #('210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_1000_after_noApod_ex','aligned','tau1000, after'),
        #('210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_3500_before_noApod_ex','prealignment','tau3500, before'),
        #('210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_3500_after_noApod_ex','aligned','tau3500, after'),
        #('210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_11_before_noApod_ex','prealignment','tau11135, before'),
        ('210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_11_after_noApod_ex','aligned','tau11135, after'),
        ]:
        d = find_file(filename, exp_type='ODNP_NMR_comp/processed',
                expno=expno)
        ph1_val = 1
        ph2_val = 0
        signal_pathway = {"ph1": ph1_val, "ph2": ph2_val}
        excluded_pathways = [(0, 0),(0,-1)]
        def select_pathway(s,pathway):
            retval = s
            for k,v in pathway.items():
                retval = retval[k,v]
            return retval    
        d.ift('t2')
        # Zeroth order correction
        ph0 = select_pathway(d['t2':0],signal_pathway)
        if len(ph0.dimlabels) > 0:
            assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
            ph0 = zeroth_order_ph(ph0)
            print('phasing dimension as one')
        else:
            print("there is only one dimension left -- standard 1D zeroth order phasing")
            ph0 = ph0/abs(ph0)
        d /= ph0
        fl.basename = expno
        # Put what goes into here in the time domain
        freq_lim = integrate_limits(select_pathway(d,signal_pathway),
                convolve_method='Lorentzian',
                fl=fl)
fl.show();quit()
