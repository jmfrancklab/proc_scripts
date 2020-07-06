from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from sympy import symbols
fl = fl_mod()
t2 = symbols('t2')


for searchstr, exp_type, nodename, postproc, label_str, slice_f in [
        ('200302_alex_probe_water','test_equip','signal','spincore_Hahn_echoph_v1','microwaves off',(-5e3,5e3)),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    #nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    print(ndshape(s))
    s.reorder('t',first=True)
    t2_axis = s.getaxis('t')[0:2048]
    s.setaxis('t',None)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.setaxis('t2',t2_axis)
    s.setaxis('nScans',r_[0:nScans])
    s.reorder('t2',first=False)
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
    logger.info(strm(ndshape(s)))
    #}}}
    #{{{ apply phase corrections
    best_shift = hermitian_function_test(s)
    logger.info(strm("best shift is",best_shift))
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
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
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    fl.next('time domain')
    fl.plot(s_sliced)
    s_sliced.ft('t2')
    fl.next('Spectrum FT')
    fl.plot(s_sliced.real, alpha=0.5, label='real - %s'%label_str)
    fl.plot(s_sliced.imag, alpha=0.5, label='imag - %s'%label_str)
fl.show();quit()
