from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
rcParams["savefig.transparent"] = True
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 5e3
color_choice = True
for date,id_string,label_str,color_str in [
        ('200302','alex_probe_water','newest','blue'),
        ('191212','echo_5','older','red'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.reorder('t2',first=False)
    fl.next('raw data -- coherence channels')
    s.ft(['ph2','ph1'])
    fl.image(s)
    fl.next('filtered + rough centered data %s'%label_str)
    s.ft('t2', shift=True)
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.image(s)
    s.ft('t2')
    k = s.C
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1],shift_val=1)
    fl.next('hermitian test')
    fl.plot(residual)
    print("best shift is",best_shift)
    # {{{ slice out the FID appropriately and phase correct
    # it
    s.mean('nScans')
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
    s = s['t2':(0,None)]
    fl.next('phased - time')
    fl.plot(s['ph2',-2]['ph1',1])
    s.ft('t2')
    #s.convolve('t2',7)
    if label_str == 'microwaves on':
        s *= -1
    fl.next('')
    s.name('')
    if color_choice:
        fl.plot(s['ph2',-2]['ph1',1],label='%s'%label_str,c=color_str)
    else:
        fl.plot(s['ph2',-2]['ph1',1],label='%s'%label_str)
fl.show();quit()
