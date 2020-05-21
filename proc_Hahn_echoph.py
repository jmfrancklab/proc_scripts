from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from sympy import symbols
fl = fl_mod()
t2 = symbols('t2')
filter_bandwidth = 5e3
slice_f = (-5e3,5e3)
color_choice = True
for searchstr,exp_type,nodename,postproc,label_str,color_str in [
        ('200302_alex_probe_water','test_equip','signal','Hahn_echoph','microwaves off','blue'),
        ]:
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict)
    s.mean('nScans')    
    fl.next('signal')
    s = s['t2':slice_f]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    logger.info(strm(ndshape(s)))
    #s = s['ph1',1]['ph2',0]
    best_shift = hermitian_function_test(s)
    logger.info(strm("best shift is",best_shift))
    # {{{ slice out the FID appropriately and phase correct
    # it
    s.ft('t2')
    s_uncorrected = s.C
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    # apply a lorentzian broadening
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    fl.next('time domain after hermitian test')
    fl.plot(s)
    ph0 = s['t2':0]['ph2',0]['ph1',1]
    logger.info(strm(ndshape(ph0)))
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0, fl=fl)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2', pad=512) # power of 2 FT
    s.convolve('t2',10) # so that resolution of plot isn't higher than that of screen
    fl.image(s)
    s = s['ph1',1]['ph2',0].C
    s.ift('t2')
    s = s['t2':(0,None)]
    s.ft('t2')
    fl.next('processed data')
    fl.plot(s_uncorrected['ph2',-2]['ph1',1],
            label='without time-axis correction',c='k')
    fl.plot(s,label='with time-axis correction',c='r')
fl.show();quit()
