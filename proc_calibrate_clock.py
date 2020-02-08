from pyspecdata import *
from align_slice import align_and_slice
fl = figlist_var()
for date,id_string in [
        ('200206','calibrate_clock_1'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    s.rename('t','t2').set_units('t2','s')
    clock_correction = 0
    centerpoint = abs(s).mean('vd').argmax('t2').item()
    s.setaxis('t2', lambda x: x-centerpoint)
    s.ft('t2', shift=True)
    s = s['t2':(-1e3,1e3)]
    fl.next('image raw')
    fl.image(s)
    fl.next('image shifted and sliced')
    s = align_and_slice(s, fl=fl)
    fl.image(s)
    fl.next('phase error vs vd')
    fl.plot(s.sum('t2').angle, 'o')
    fl.next('phase error, unwrapped vs vd')
    s = s['vd',1:]/s['vd',:-1]
    s = s.angle.name('signal phase').set_units('rad')
    s.data = s.data.cumsum()
    fl.plot(s,'o')
    #print (s.data[-1]-s.data[-2])/(s.getaxis('vd')[-1]-s.getaxis('vd')[-2])
fl.show()
