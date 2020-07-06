from pyspecdata import *
from sympy import symbols
from proc_scripts import * 
def center_CPMG_echo(s):
    s.register_axis({'t2':0})
    s /= zeroth_order_ph(s['t2':0])
    time_bound = min(abs(s.getaxis('t2')[r_[0,-1]]))
    s = s['t2':(-time_bound,time_bound)]
    return s
    
