from pyspecdata import *
from scipy.optimize import minimize
from proc_scripts import hermitian_function_test,zeroth_order_ph,postproc_dict,fl_mod
from sympy import symbols
from numpy import *
fl = fl_mod()
t2 = symbols('t2')
logger = init_logging("info")
for searchstr,exp_type,nodename,postproc in [
    ['201211_Ni_sol_probe_nutation_1','nutation','nutation','spincore_nutation_v1']
    ]:
    s = find_file(searchstr,exp_type=exp_type,expno=nodename,postproc=postproc,
            lookup=postproc_dict)#,fl=fl) 
    #s = s['t2':(-1000,1000)]
    #fl.show();quit()
    # {{{ do the rough centering before anything else!
    # in particular -- if you don't do this before convolution, the
    # convolution doesn't work properly!
    s.ift('t2')
    # }}}
    s.setaxis('t2', lambda x: x-abs(s['ph1',1]['ph2',-2]).mean_all_but('t2').argmax('t2').item())
    fl.next('t domain after removing 1st couple pts')
    #s = s['t2',2:]
    fl.image(s)
    fl.next('freq domain after filter')
    s.ft('t2')
    fl.image(s)
    #fl.show();quit()
    s.ift('t2')
    # {{{ centering of data using hermitian function test
    best_shift = hermitian_function_test(s['ph2',0]['ph1',1])
    logger.info(strm("best shift is",best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('with time shift')
    fl.image(s)
    fl.next('freq domain after time correction')
    s.ft('t2')
    fl.image(s)
    #fl.show();quit()
    #}}}
    
    #s.ft('t2',pad=4096)
    
    #{{{ selecting coherence and convolving
    s = s['ph2',0]['ph1',1]
    fl.next('select $\\Delta p$ and convolve')
    s.convolve('t2',50)
    fl.image(s)
    #fl.show();quit()
    #}}}
    
    #{{{ slicing
    s = s['t2':(-5000,13000)]
    fl.next('sliced')
    fl.image(s)
    #fl.show();quit()
    #}}}
    
    #{{{ phasing with zeroth order correction
    s.ift('t2')
    fl.next('final time domain')
    ph0 = zeroth_order_ph(s['t2':0], fl=None)
    s /= ph0
    fl.image(s)
    fl.next('phased')
    s.ft('t2',pad=4096)
    fl.image(s)
    fl.real_imag('phased data',s)
    #}}}
fl.show();quit()
