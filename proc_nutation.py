from pyspecdata import *
from scipy.optimize import minimize
from hermitian_function_test import zeroth_order_ph
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
date = '200219'
for id_string in [
    'nutation_alex_probe',
    ]:
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
    s.ft('t2',shift=True,pad=4096)
    #s *= exp(1j*2*pi*0.42) # manually determined ph correction
    fl.next('image, raw')
    fl.image(s)
    s.ft(['ph2','ph1'])
    fl.next('image, all coherence channels')
    fl.image(s)
    s = s['ph2',0]['ph1',1].C
    fl.next(id_string+'image')
    fl.image(s)
    s.ift('t2',pad=4096)
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    s.ft('t2',pad=4096)
    s = s['t2':(-100,50)]
    #s.ft('t2',pad=4096)
    fl.next('sliced')
    fl.image(s)
    s.ift('t2')
    ph0 = zeroth_order_ph(s, fl=None)
    s /= ph0
    fl.next('phased')
    s.ft('t2',pad=4096)
    fl.image(s)
    fl.show();quit()
            
    fl.next(id_string+'image -- $B_1$ distribution')
    fl.image(abs(s.C.ft('p_90',shift=True)))
    fl.next(id_string+'image abs')
    fl.image(abs(s))
fl.show();quit()
