from hermitian_function_test import zeroth_order_ph,hermitian_function_test
from proc_scripts import * 
from pyspecdata import *
from sympy import symbols

def slice_FID_from_echo(s,time_range):
    best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1], fl=fl)
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    ph0 = s['t2':0]['ph2',-2]['ph1',1]
    s /= zeroth_order_ph(ph0, fl=None)
    if s['t2':0]['ph2',-2]['ph1',1]['power',0].real < 0:
        s *= -1
    fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range',
        s,time_range)
    s = s['t2':time_range]
    fl.next('echo mirror test')
    echo_start = s.getaxis('t2')[0]
    dw = diff(s.getaxis('t2')[r_[0,1]]).item()
    centered_echo = s['t2':(echo_start,-echo_start+dw)]
    plotdata = abs(centered_echo/centered_echo['t2',::-1])
    plotdata[lambda x: x>2] = 2
    fl.image(plotdata)
    fl.next('apodize and zero fill')
    R = 5.0/(time_range[-1]) # assume full decay by end time
    s *= exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1024)
    fl.image(s)
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    return s
