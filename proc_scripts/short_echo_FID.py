from hermitian_function_test import zeroth_order_ph
from proc_scripts import side_by_side
from pyspecdata import *
def FID(s):
    ph0 = s['t2':0]['ph2',-2]['ph1',1]
    s /= zeroth_order_ph(ph0, fl=fl)
    if s['t2':0]['ph2',-2]['ph1',1]['power',0].real < 0:
        s *= -1
    fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range',
            s,time_range)
    s = s['t2':time_range]
    # {{{ all of this is to check and see if we think
    # we can add the two sides of the echo to increase
    # our SNR -- right now, it doesn't look like it
    fl.next('echo mirror test')
    echo_start = s.getaxis('t2')[0]
    dw = diff(s.getaxis('t2')[r_[0,1]]).item()
    centered_echo = s['t2':(echo_start,-echo_start+dw)]
    plotdata = abs(centered_echo/centered_echo['t2',::-1])
    plotdata[lambda x: x>2] = 2
    fl.image(plotdata)
    # }}}
    fl.next('apodize and zero fill')
    R = 5.0/(time_range[-1]) # assume full decay by end time
    s *= exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1024)
    fl.image(s)
    # {{{ select FID
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    return s
