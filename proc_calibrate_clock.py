from pyspecdata import *
from align_slice import align_and_slice
from scipy.optimize import leastsq
fl = figlist_var()
apply_correction = True
for date,id_string in [
        ('200210','calibrate_clock_1'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    s.rename('t','t2').set_units('t2','s')
    s.setaxis('vd',s.getaxis('vd')*1e-6) # not sure why, but vd list needed to be saved in us instead of sec
    centerpoint = abs(s).mean('nScans').mean('vd').argmax('t2').item()
    s.setaxis('t2', lambda x: x-centerpoint)
    fl.next('image raw -- time domain')
    fl.image(s, interpolation='bicubic')
    s.ft('t2', shift=True)
    s = s['t2':(-5e3,5e3)]
    if apply_correction:
        clock_correction = 1.5486
        s *= exp(-1j*s.fromaxis('vd')*clock_correction)
    fl.next('image raw')
    fl.image(s)
    fl.next('image shifted and sliced')
    s = align_and_slice(s,convwidth=0,fl=fl)
    s.mean('nScans')
    fl.image(s)
    fl.next('phase error vs vd')
    fl.plot(s.sum('t2').angle, 'o')
    fl.next('phase error, unwrapped vs vd - sig avg')
    s = s['vd',1:]/s['vd',:-1]
    s = s.angle.name('signal phase').set_units('rad')
    s.data = s.data.cumsum()
    fl.plot(s,'o',human_units=False)
    # begin fit to return clock correction
    x = s.getaxis('vd')
    fitfunc = lambda p, x: p[0]*x
    errfunc = lambda p_arg,x_arg,y_arg: fitfunc(p_arg, x_arg) - y_arg
    p_ini = [1.0]
    p1,success = leastsq(errfunc,p_ini[:],args=(x,s.data))
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit,fitfunc(p1,x_fit),':',human_units=False,label='%f'%p1[0])
    clock_correction = p1[0]
    print(clock_correction)
fl.show()
