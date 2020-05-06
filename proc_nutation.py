from pyspecdata import *
from scipy.optimize import minimize
from proc_scripts import hermitian_function_test,zeroth_order_ph,load_data
from sympy import symbols
from numpy import *
fl = figlist_var()
t2 = symbols('t2')

for searchstr,exp_type,nodename in [
    ['200219_nutation_alex_probe','test_equip','nutation']
    ]:
    s = load_data(searchstr,exp_type=exp_type,which_exp=nodename,postproc='nutation') 
    # {{{ do the rough centering before anything else!
    # in particular -- if you don't do this before convolution, the
    # convolution doesn't work properly!
    s.ft(['ph2','ph1'])
    rough_center = abs(s)['ph2',0]['ph1',1].convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.next('raw')
    fl.image(s)
    print(ndshape(s))
    # }}}
    # {{{ centering of data using hermitian function test
    best_shift = hermitian_function_test(s['ph2',0]['ph1',1])
    fl.next('hermitian test')
    fl.plot(best_shift)
    print("best shift is",best_shift)
    s.ft('t2', shift=True)
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    fl.next('time domain after hermitian test')
    #fl.image(s)
    #fl.show();quit()
    #}}}
    #{{{ reviewing data imaging thus far
    fl.next('time domain (all $\\Delta p$)')
    fl.image(s)
    fl.next('frequency domain (all $\\Delta p$)')
    s.ft('t2',pad=4096)
    fl.image(s)
    #fl.show();quit()
    #}}}
    #{{{ selecting coherence and convolving
    s = s['ph2',0]['ph1',1].C
    fl.next('select $\\Delta p$ and convolve')
    s.convolve('t2',50)
    fl.image(s)
    fl.next('Figure 1')
    fl.image(s)
    #fl.show();quit()
    #}}}
    #{{{ slicing
    s = s['t2':(-250,250)]
    #s.ft('t2',pad=4096)
    fl.next('sliced')
    fl.image(s)
    #}}}
    #{{{ phasing with zeroth order correction
    s.ift('t2')
    fl.next('final time domain')
    ph0 = zeroth_order_ph(s['t2':0], fl=None)
    s /= ph0
    fl.image(s)
    #fl.show();quit()
    fl.next('phased')
    s.ft('t2',pad=4096)
    fl.image(s)
    fl.next('real')
    fl.image(s.real)
    my_clim = gci().get_clim()
    fl.next('imag')
    fl.image(s.imag)
    gci().set_clim(my_clim) # to match real
    #}}}
    #        
    #fl.next(id_string+'image -- $B_1$ distribution')
    #fl.image(abs(s.C.ft('p_90',shift=True)))
    #fl.next(id_string+'image abs')
    #fl.image(abs(s))
fl.show();quit()
