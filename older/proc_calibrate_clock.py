from pyspecdata import *
from align_slice import align_and_slice
fl = figlist_var()
apply_correction = True
for date,id_string in [
        ('200212','calibrate_clock_6'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    s.rename('t','t2').set_units('t2','s')
    centerpoint = abs(s).mean('nScans').mean('vd').argmax('t2').item()
    s.setaxis('t2', lambda x: x-centerpoint)
    fl.next('image raw -- time domain')
    fl.image(s, interpolation='bicubic')
    s.ft('t2', shift=True)
    s = s['t2':(-0.15e3,0.15e3)]
    if apply_correction:
        #clock_correction = 1.692
        #clock_correction = -1.089
        clock_correction = 1.785
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
    raise RuntimeError("the following code should not be using leastsq -- this is a line, while leastsq is for non-linear fits")
    fitfunc = lambda p, x: p[0]*x
    errfunc = lambda p_arg,x_arg,y_arg: fitfunc(p_arg, x_arg) - y_arg
    p_ini = [1.0]
    p1,success = leastsq(errfunc,p_ini[:],args=(x[-4:],s.data[-4:]))
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit,fitfunc(p1,x_fit),':',human_units=False,label='%f'%p1[0])
    clock_correction = p1[0]
    print(clock_correction)
fl.show()
