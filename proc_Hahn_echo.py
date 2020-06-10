from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
for date,id_string in [
        ('200219','alex_echo_AGprobe_1'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip'))
    print(ndshape(s))
    nEchoes = s.get_prop('acq_params')['nEchoes']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    s.set_units('t','s')
    print((s.get_prop('acq_params')))
    print((s.get_prop('nScans')))
    fl.next('raw data')
    fl.plot(s.real,alpha=0.4)
    #fl.plot(s.imag,alpha=0.4)
    #fl.plot(abs(s),':',c='k',alpha=0.4)
    s.ft('t',shift=True)
    s.convolve('t',10)
    fl.next('Absolute value -- FT')
    #fl.plot(s.real,alpha=0.4)
    #fl.plot(s.imag,alpha=0.4)
    fl.plot(abs(s),c='red')
    fl.next('Real - FT')
    fl.plot(s.real,alpha=0.4)
    #fl.plot(s.imag,alpha=0.4)
    #fl.plot(abs(s),':',c='k',alpha=0.4)
fl.show()
