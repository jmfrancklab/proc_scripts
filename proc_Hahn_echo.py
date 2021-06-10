from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
phase_cycled=True
for date,id_string in [
        ('210318','TEMPOL500mM_cap_probe_echo_2'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip'))
    nEchoes = s.get_prop('acq_params')['nEchoes']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    s.set_units('t','s')
    if phase_cycled:
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
        s.setaxis('ph2',r_[0.,2.]/4)
        s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
        fl.next('imag')
        s.mean('nScans')
        fl.image(s)
        s.ft('t2',shift=True)
        fl.next('image-FT')
        fl.image(s)
        fl.next('coherence,FT')
        s.ft(['ph2','ph1'])
        fl.image(s)
        fl.next('data plot')
        fl.plot(s['ph1',1]['ph2',0])
        fl.plot(s.imag['ph1',1]['ph2',0])
    else:
        s.ft('t2',shift=True)
        s.ft(['ph1','ph2'])
    fl.next('coherence domain')
    fl.image(s)
    s.mean_all_but('t2')
    fl.next('Absolute value -- FT for capillary probe')
    fl.plot(s.real,alpha=0.4,label='real')
    fl.plot(s.imag,alpha=0.4,label='imaginary')
    fl.next('Real - FT')
    fl.plot(s.real,alpha=0.4,label='real')
    fl.plot(s.imag,alpha=0.4,label='imaginary')
fl.show()
