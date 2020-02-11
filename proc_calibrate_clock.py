from pyspecdata import *
from align_slice import align_and_slice
from scipy.optimize import leastsq
fl = figlist_var()
for date,id_string in [
        ('200210','calibrate_clock_1'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    s.rename('t','t2').set_units('t2','s')
    clock_correction = 0
    #clock_correction = (1.46491/999104e-6)
    clock_correction = 1.54862309
    centerpoint = abs(s).mean('nScans').mean('vd').argmax('t2').item()
    s.setaxis('t2', lambda x: x-centerpoint)
    fl.next('image raw -- time domain')
    fl.image(s, interpolation='bicubic')
    s.ft('t2', shift=True)
    s = s['t2':(-5e3,5e3)]
    s *= exp(-1j*s.fromaxis('vd')*1e-6*clock_correction)
    fl.next('image raw')
    fl.image(s)
    fl.next('image shifted and sliced')
    s = align_and_slice(s,convwidth=0,fl=fl)
    s.mean('nScans')
    fl.image(s)
    fl.show();quit()
    fl.next('phase error vs vd')
    fl.plot(s.sum('t2').angle, 'o')
    fl.next('phase error, unwrapped vs vd - sig avg')
    s = s['vd',1:]/s['vd',:-1]
    s = s.angle.name('signal phase').set_units('rad')
    s.data = s.data.cumsum()
    s.setaxis('vd',s.getaxis('vd')*1e-6)
    fl.plot(s,'o',human_units=False)
    # begin fit to return clock correction
    x = s.getaxis('vd')
    fitfunc = lambda p, x: p[0]*x
    errfunc = lambda p_arg,x_arg,y_arg: fitfunc(p_arg, x_arg) - y_arg
    p_ini = [1.0]
    p1,success = leastsq(errfunc,p_ini[:],args=(x,s.data))
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit,fitfunc(p1,x_fit),':',human_units=False)
    print(p1)
fl.show()
