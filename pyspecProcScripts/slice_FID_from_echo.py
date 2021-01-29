from proc_scripts import * 
from pyspecdata import *
from sympy import symbols

def slice_FID_from_echo(s, ph1=1, ph2=-2):
    fl = fl_mod() 
    best_shift = hermitian_function_test(s[
        'ph2',ph2]['ph1',ph1])
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    coh_slice = s['t2':0]['ph2',ph2]['ph1',ph1]
    if len(coh_slice.dimlabels) > 0:
        assert len(coh_slice.dimlabels) == 1, repr(ndshape(coh_slice.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(coh_slice, fl=fl)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = coh_slice/abs(coh_slice)
    s /= ph0
    s['t2',0] *= 0.5
    if 'power' in s.dimlabels:
        if s['t2':0]['ph2',ph2]['ph1',ph1]['power',0].real < 0:
            s *= -1
    elif 'vd' in s.dimlabels:
        if s['t2':0]['ph2',ph2]['ph1',ph1]['vd',-1].real < 0:
            s *= -1
    return s
