from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import load_data
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 5e3
color_choice = True
for searchstr,label_str,color_str in [
        ('191017_echo_adc44','microwaves off','blue'),
        ]:
    #nPhaseSteps = 8
    s = load_data(searchstr,'Hahn_echoph',None)
    s.ft('t2',shift=True)
    fl.next('raw data, chunked')
    fl.image(abs(s))
    s.ft(['ph1','ph2'])
    fl.next('coherence')
    fl.image(abs(s))
    s = s['ph1',1]['ph2',0].C
    s.mean('nScans')#,return_error=False)
    fl.next('signal')
    fl.plot(abs(s),label=label_str)
    slice_f = (-5e3,5e3)
    s = s['t2':slice_f].C
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    s.mean('nScans')
    s.ft('t2')
    fl.image(s)
    #k = s.C*exp(-1j*s.fromaxis('t2')*0.9*2*pi)
    k = s.C*exp(-1j*0.9*2*pi)
    k *= exp(-1j*k.fromaxis('t2')*2*pi*0.005)
    s = k.C
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1],shift_val=1)
    fl.next('hermitian test')
    fl.plot(residual)
    print("best shift is",best_shift)
    # {{{ slice out the FID appropriately and phase correct
    # it
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    s *= exp(-s.getaxis('t2')/40e-3)
    fl.next('time domain after hermitian test')
    fl.image(s)
    ph0 = s['t2':0]['ph2',-2]['ph1',1]
    print(ndshape(ph0))
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0, fl=fl)
        print('phasing dimension as one')
    else:
        print("there is only one dimension left -- standard 1D zeroth order phasing")
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    k.ift('t2')
    s = s['t2':(0,None)]
    k = k['t2':(0,None)]
    fl.next('phased - time')
    fl.plot(s['ph2',-2]['ph1',1])
    s.ft('t2')
    k.ft('t2')
    #s.convolve('t2',7)
    fl.next('')
    s.name('')
    k.rename('t2','Offset').set_units('Offset','Hz')
    s.rename('t2','Offset').set_units('Offset','Hz')
    fl.plot(k['ph2',-2]['ph1',1],label='without time-axis correction',c='k')
    fl.plot(s['ph2',-2]['ph1',1],label='with time-axis correction',c='r')
fl.show();quit()
