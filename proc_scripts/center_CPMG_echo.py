from pyspecdata import *
from sympy import symbols
from proc_scripts import * 
def center_CPMG_echo(s):
    best_shift = hermitian_function_test(s)
    logger.info(strm("best shift is",best_shift))
    s.setaxis('t2', lambda x: x-best_shift)
    s.register_axis({'t2':0})
    s /= zeroth_order_ph(s['t2':0])
    time_bound = min(abs(s.getaxis('t2')[r_[0,-1]]))
    s = s['t2':(-time_bound,time_bound)]
    assert isclose(s.getaxis('t2')[0],-s.getaxis('t2')[-1]),"echo is not symmetric! you are using the wrong code!!"
    return s
