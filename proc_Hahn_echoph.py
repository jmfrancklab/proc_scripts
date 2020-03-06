from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 5e3
for date,id_string,label_str in [
        ('191007','echo_2','n'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = 8 
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    print(ndshape(s))
    s.reorder('t',first=True)
    t2_axis = s.getaxis('t')[0:2048]
    s.setaxis('t',None)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.set_units('t2','s')
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    fl.next('ft')
    fl.plot(s)
    s = s['t2':(-2e3,2e3)]
    s.ift('t2')
    fl.next('raw data')
    fl.image(s)
    fl.next('raw data -- coherence pathways')
    #s = s['power',0].C
    s.ft(['ph2','ph1'])
    #s *= exp(-s.getaxis('t2')/10e-3)
    fl.image(s)
    fl.show();quit()
    fl.next('filtered + rough centered data')
    s.ft('t2', shift=True)
    s *= exp(-1j*2*pi*10.3)
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.image(s)
    s.ft('t2')
    k = s.C
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1])
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
    s.ft('t2',pad=512)
    fl.image(s)
    s.ift('t2')
    s = s['t2':(0,None)]
    fl.next('phased - time')
    fl.plot(s)
    s.ft('t2')
    fl.next('phased')
    s.name('')
    fl.plot(k['ph2',-2]['ph1',1],c='k')
    fl.plot(s['ph2',-2]['ph1',1])
    fl.plot(s.imag['ph2',-2]['ph1',1])
    fl.show();quit()
