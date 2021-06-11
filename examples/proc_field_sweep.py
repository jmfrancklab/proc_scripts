"""
Field Sweep
===========

A field sweep to detect the ESR via NMR-ODNP!
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from scipy.optimize import leastsq,minimize,basinhopping
from sympy import symbols
import pywt
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 20e3
filename = '210611_S175R1a_pR_DDM_field_dep'
for nodename,postproc,label_str,freq_slice,field_slice in [
        ('32dBm_finer','field_sweep','Sams field sweep',(-700,700),(-400,300)),
        ]:
    s = find_file(filename,exp_type='ODNP_NMR_comp/field_dependent',
            expno=nodename,postproc=postproc,lookup=postproc_dict,fl=fl)
    s = s['t2':freq_slice]
    if s.get_prop('acq_params')['nPhaseSteps'] == 8:
        s.mean('nScans')
        s = s['ph1',1]['ph0',0].C
    else:
        s=s['ph1',1]['power',0]['nScans',0].C
    fl.next('filtered + rough centered data')
    s.ift('t2')
    s.ft('t2')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).C.convolve('t2',0.0001).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    s.ft('t2')
    fl.next('line plots')
    for z in range(len(s.getaxis('Field'))):
        fl.plot(abs(s['Field',z]),label='%d'%z)
    s_ = s['t2':field_slice].sum('t2')
    fl.next('sweep, without hermitian')
    fl.plot(abs(s_),'o-')
    fl.show();quit()    
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s['ph2',-2]['ph1',1])
    fl.next('hermitian test')
    fl.plot(residual)
    print("best shift is",best_shift)
    # {{{ slice out the FID appropriately and phase correct
    # it
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    fl.next('time domain after hermitian test')
    fl.image(s)
    for field_index,this_field in enumerate(s.getaxis('Field')):
        ph0 = s['t2':0]['ph2',-2]['ph1',1]['Field',field_index]
        if len(ph0.dimlabels) > 1:
            assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
            ph0 = zeroth_order_ph(ph0, fl=fl)
            print('phasing dimension as one')
        else:
            print("there is only one dimension left -- standard 1D zeroth order phasing")
            ph0 = ph0/abs(ph0)
        s['Field',field_index] /= ph0
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.image(s)
    #fl.show();quit()
    s = s['ph2',-2]['ph1',1]
    s = s['t2':(-1e3,1e3)]
    print(ndshape(s))
    fl.next('AOT RM 4AT, capillary probe\nPhased (individual loaded files) - 0.833 mM net, 50 mM local')
    fl.plot(s,alpha=0.6,label='%s'%label_str)
    fl.plot(s.imag,alpha=0.6,label='%s'%label_str)
    s = s['t2':(-200,400)].C.sum('t2')
    fl.next('Field dependent')
    fl.plot(s,'o')
    print(s.getaxis('Field'))
fl.show()
