from pylab import *
from pyspecdata import *
from scipy.optimize import minimize
from proc_scripts import hermitian_function_test,zeroth_order_ph,postproc_dict,fl_mod
from sympy import symbols
from numpy import *
fl = fl_mod()
t2 = symbols('t2')
logger = init_logging("info")
max_kHz = 200
for searchstr,exp_type,nodename,postproc,freq_slice in [
    #['201211_Ni_sol_probe_nutation_1','nutation','nutation',
    #    'spincore_nutation_v1',(-5000,13000)],
    #['210302_Ni_cap_probe_nutation_1','nutation','nutation',
    #    'spincore_nutation_v1',(-4e3,1.2e4)]
    ['210120_Ni_sol_probe_nutation_amp_3','nutation','nutation',
        'spincore_nutation_v2',(-15.2e3,11.01e3)]
    ]:
    s = find_file(searchstr,exp_type=exp_type,expno=nodename,postproc=postproc,
            lookup=postproc_dict) 
    if 'amp' in s.dimlabels:
        plen = s.get_prop('acq_params')['p90_us']*1e-6
        logger.info(strm('pulse length is:',plen))
    # {{{ do the rough centering before anything else!
    # in particular -- if you don't do this before convolution, the
    # convolution doesn't work properly!
    s.ift('t2')
    # {{{ centering of data using hermitian function test
    else:
        best_shift = hermitian_function_test(s['ph2',0]['ph1',1])
        logger.info(strm("best shift is",best_shift))
        s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
        fl.next('with time shift')
        fl.image(s)
        fl.next('freq domain after time correction')
        s.ft('t2')
        fl.image(s)
    #}}}
    if 'amp' in s.dimlabels:
        s.setaxis('t2', lambda x: x-abs(s['ph1',1]['ph2',-2]).mean_all_but('t2').argmax('t2').item())
    #}}}
    fl.next('t domain centered')
    fl.image(s)
    fl.next('freq domain centered')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    #{{{ selecting coherence and convolving
    s = s['ph2',0]['ph1',1]
    fl.next('select $\\Delta p$')
    else:
        s.convolve('t2',50)
    fl.image(s)
    #}}}
    #{{{ slicing
    s = s['t2':freq_slice]
    fl.next('sliced')
    fl.image(s)
    #}}}
    if 'amp' in s.dimlabels:
        s.setaxis('amp',lambda x:x*plen)
        s.set_units('amp','s')
        ind_dim = '\\tau_p a'
        s.rename('amp',ind_dim)
    else:
        ind_dim = 'p_90'
    #{{{ phasing with zeroth order correction
    fl.next('final time domain')
    ph0 = zeroth_order_ph(s['t2':0], fl=None)
    s /= ph0
    fl.image(s)
    fl.next('phased')
    fl.image(s)
    fl.real_imag('phased data',s)
    #}}}
    fl.next('FT')
    title('FT to get $\gamma B_1/a$')
    s.ft(ind_dim,shift=True)
    fl.image(s[ind_dim:(-1e3*max_kHz,1e3*max_kHz)])
    fl.next('absFT')
    title('FT to get $\gamma B_1/a$')
    fl.image(abs(s[ind_dim:(-1e3*max_kHz,1e3*max_kHz)]))
    gridandtick(gca(),gridcolor=[1,1,1])
fl.show();quit()

