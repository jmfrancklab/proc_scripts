from pyspecdata import *
fl = figlist_var()
date = '200206'
for id_string in [
        ('calibrate_clock_1'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    s.rename('t','t2').set_units('t2','s')
    clock_correction = 0
    fl.next('image raw')
    acq_time = s.getaxis('t2')[-1]
    manual_taxis_zero = acq_time/2.0
    s.setaxis('t2',lambda x: x-manual_taxis_zero)
    s.ft('t2', shift=True)
    data = s.C
    fl.image(s)
    fl.next('phase error vs vd')
    fl.plot(s.sum('t2').angle, 'o')
    fl.next('phase error, unwrapped vs vd')
    s = s['vd',1:]/s['vd',:-1]
    s = s.angle.name('signal phase').set_units('rad')
    s.data = s.C.data.cumsum()
    fl.plot(s,'o')
    print((s.data[-1]-s.data[-2])/(s.getaxis('vd')[-1]-s.getaxis('vd')[-2]))
    clock_correction = 9.35
    fl.next('image corrected')
    data *= exp(1j*data.fromaxis('vd')*clock_correction)
    fl.image(data)
fl.show()
