from pyspecdata import *
from scipy.optimize import minimize
from hermitian_function_test import hermitian_function_test,zeroth_order_ph
from sympy import symbols
from numpy import *
fl = figlist_var()
t2 = symbols('t2')

date = '200219'
for id_string in [
    'nutation_alex_probe',
    ]:
    #{{{ Loading in data
    filename = date+'_'+id_string+'.h5'
    nodename = 'nutation'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    orig_t = s.getaxis('t')
    s.set_units('p_90','s')
    #s.setaxis('t',None)
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.reorder('t2',first=False)
    #}}}
    # {{{ do the rough centering before anything else!
    # in particular -- if you don't do this before convolution, the
    # convolution doesn't work properly!
    s.ft(['ph2','ph1'])
    rough_center = abs(s)['ph2',0]['ph1',1].convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
        
    # }}}
    # {{{ centering of data using hermitian function test
    residual,best_shift = hermitian_function_test(s['ph2',0]['ph1',1])
    fl.next('hermitian test')
    fl.plot(residual)
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
    
    #fl.next('imag')
    #fl.image(s.imag)
    #gci().set_clim(my_clim) # to match real
    #}}}
    #        
    #fl.next(id_string+'image -- $B_1$ distribution')
    #fl.image(abs(s.C.ft('p_90',shift=True)))
    #fl.next(id_string+'image abs')
    #fl.image(abs(s))
fl.show();quit()
