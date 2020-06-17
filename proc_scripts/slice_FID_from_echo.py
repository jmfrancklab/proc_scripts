from proc_scripts import * 
from pyspecdata import *
from sympy import symbols

def slice_FID_from_echo(s,ph2,ph1):
    fl = fl_mod() 
    best_shift = hermitian_function_test(s[
        'ph2',ph2]['ph1',ph1])
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    ph0 = s['t2':0]['ph2',ph2]['ph1',ph1]
    s /= zeroth_order_ph(ph0, fl=None)
    if 'power' in s.dimlabels:
        if s['t2':0]['ph2',ph2]['ph1',ph1]['power',0].real < 0:
            s *= -1
    return s
