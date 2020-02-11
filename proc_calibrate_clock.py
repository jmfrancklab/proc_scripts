from pyspecdata import *
fl = figlist_var()
date = '200210'
#date = '200206'
for id_string in [
        ('calibrate_clock_1'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    s.rename('t','t2').set_units('t2','s')
    clock_correction = 0
    #clock_correction = -10.51/6
    clock_correction = -0.776915
    fl.next('image raw')
    fl.image(s)
    fl.next('image raw - FT')
    fl.image(s.C.ft('t2',shift=True))
    fl.next('image -- scans avg')
    fl.image(s.C.mean('nScans'))
    #fl.image(s)
    fl.next('image -- scans avg FT')
    fl.image(s.C.mean('nScans').ft('t2', shift=True))
    #fl.image(s.C.ft('t2', shift=True))
    acq_time = s.getaxis('t2')[-1]
    manual_taxis_zero = acq_time/2.0
    s.setaxis('t2',lambda x: x-manual_taxis_zero)
    s.ft('t2', shift=True)
    s *= exp(-1j*s.fromaxis('vd')*clock_correction)
    fl.next('image - sig avg')
    fl.image(s.mean('nScans'))
    #fl.image(s)
    fl.next('phase error vs vd - sig avg')
    fl.plot(s.sum('t2').angle, 'o')
    fl.next('phase error, unwrapped vs vd - sig avg')
    s = s['vd',1:]/s['vd',:-1]
    s = s.angle.name('signal phase').set_units('rad')
    s.data = s.data.cumsum()
    fl.plot(s,'o')
    print (s.data[-1]-s.data[-2])/(s.getaxis('vd')[-1]-s.getaxis('vd')[-2])
fl.show()
