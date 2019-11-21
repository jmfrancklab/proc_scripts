from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
fl = figlist_var()
for date,id_string in [
    ('191118','echo_DNP_3'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    nPoints = s.get_prop('acq_params')['nPoints']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    s.set_units('t','s')
    print s.get_prop('meter_powers')
    print ndshape(s)
    fl.next(id_string+'raw data ')
    fl.image(s)
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    t2_axis = orig_t[nPoints]
    s.setaxis('t',None)
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.setaxis('t2',t2_axis)
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    fl.next('raw data')
    fl.image(s)
    s = s['ph1',1]['ph2',0].C
    s.reorder('t2',first=True)
    s.ift('t2')
    fl.next('signal, time domain')
    fl.plot(s)
    t2_max = zeros_like(s.getaxis('power'))
    for x in xrange(len(s.getaxis('power'))):
        t2_max[x] = abs(s['power',x]).argmax('t2',raw_index=True).data
    print t2_max
    s.setaxis('t2',lambda t: t - s.getaxis('t2')[int(t2_max.mean())])
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    fl.next('signal, time domain shifted')
    fl.plot(s)
    s.ft('t2')
    #{{{ trick to try and keep sign of enhanced signal
    remember_sign = zeros_like(s.getaxis('power'))
    for x in xrange(len(s.getaxis('power'))):
        if s['power',x].data.real.sum() > 0:
            remember_sign[x] = 1.0
        else:
            remember_sign[x] = -1.0
    #}}}
    for x in xrange(len(s.getaxis('power'))):
        temp = s['power',x].C
        plot_list = [0,5,10,15]
        if x in plot_list:
            fl.next('signal, pre-phasing')
            fl.plot(temp.real, alpha=0.5, label='index %d real'%x)
            fl.plot(temp.imag, alpha=0.5, label='index %d imag'%x)
        SW = diff(temp.getaxis('t2')[r_[0,-1]]).item()
        thisph1 = nddata(r_[-4:4:2048j]/SW,'phi1').set_units('phi1','s')
        phase_test_r = temp * exp(-1j*2*pi*thisph1*temp.fromaxis('t2'))
        phase_test_rph0 = phase_test_r.C.sum('t2')
        phase_test_rph0 /= abs(phase_test_rph0)
        phase_test_r /= phase_test_rph0
        cost = abs(phase_test_r.real).sum('t2')
        ph1_opt = cost.argmin('phi1').data
        temp *= exp(-1j*2*pi*ph1_opt*temp.fromaxis('t2'))
        s['power',x] *= exp(-1j*2*pi*ph1_opt*temp.fromaxis('t2'))
        ph0 = temp.C.sum('t2')
        ph0 /= abs(ph0)
        temp /= ph0
        s['power',x] /= ph0
        if x in plot_list:
            fl.next('signal, post-phasing')
            fl.plot(temp.real, alpha=0.5, label='index %d real'%x)
            fl.plot(temp.imag, alpha=0.5, label='index %d imag'%x)
    s *= remember_sign
    fl.next('signal, phased')
    fl.plot(s)
    enhancement = s['t2':(-1e3,1e3)].C
    fl.next('Check slice')
    fl.plot(enhancement.C)
    enhancement.sum('t2').real
    enhanced = enhancement.data[1:]
    enhanced /= max(enhanced)
    fl.next('150uL TEMPOL enhancement curve')
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
    fl.plot(power_axis_W,enhanced,'.')
fl.show();show()
