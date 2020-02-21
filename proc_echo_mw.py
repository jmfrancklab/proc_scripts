from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from scipy import interpolate
fl = figlist_var()
k_sigma = True
for date,id_string in [
    #('200127','echo_DNP_TCM51C_1'),
    #('200128','echo_DNP_TCM118C_1'),
    #('200130','echo_DNP_1'),
    #('200130','echo_DNP_2'),
    #('200130','echo_DNP_3'),
    #('191118','echo_DNP_3'),
    #('191217','echo_DNP_1'),
    #('200130','echo_DNP_5'),
    #('200130','echo_DNP_AG'),
    ('200220','DNP_Y191R1apR_1'),
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
    fl.next('raw data, chunked')
    fl.image(s)
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    fl.next('coherence levels')
    fl.image(s)
    s = s['ph1',1]['ph2',0].C
    fl.next('viz - signal')
    fl.image(s)
    s.reorder('t2',first=True)
    s.ift('t2')
    t2_max = zeros_like(s.getaxis('power'))
    for x in range(len(s.getaxis('power'))):
        t2_max[x] = abs(s['power',x]).argmax('t2',raw_index=True).data
    s.setaxis('t2',lambda t: t - s.getaxis('t2')[int(t2_max.mean())])
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    fl.next('t=0, FID, then FT')
    fl.plot(s)
    fl.next('viz - signal 2')
    fl.image(s)
    remember_sign = zeros_like(s.getaxis('power'))
    for x in range(len(s.getaxis('power'))):
        if s['power',x].data.real.sum() > 0:
            remember_sign[x] = 1.0
        else:
            remember_sign[x] = -1.0
    temp = s['power',-4].C
    fl.next('signal, comparison')
    fl.plot(temp.real, alpha=0.5, label='real, pre-phasing')
    fl.plot(temp.imag, alpha=0.5, label='imag, pre-phasing')
    SW = diff(temp.getaxis('t2')[r_[0,-1]]).item()
    thisph1 = nddata(r_[-4:4:5000j]/SW,'phi1').set_units('phi1','s')
    phase_test_r = temp * exp(-1j*2*pi*thisph1*temp.fromaxis('t2'))
    phase_test_rph0 = phase_test_r.C.sum('t2')
    phase_test_rph0 /= abs(phase_test_rph0)
    phase_test_r /= phase_test_rph0
    cost = abs(phase_test_r.real).sum('t2')
    ph1_opt = cost.argmin('phi1').data
    temp *= exp(-1j*2*pi*ph1_opt*temp.fromaxis('t2'))
    s *= exp(-1j*2*pi*ph1_opt*temp.fromaxis('t2'))
    ph0 = temp.C.sum('t2')
    ph0 /= abs(ph0)
    temp /= ph0
    s /= ph0
    fl.next('signal, comparison')
    fl.plot(temp.real, alpha=0.5, label='real, post-phasing')
    fl.plot(temp.imag, alpha=0.5, label='imag, post-phasing')
    # for some reason, signs are exactly inverted when phased this way
    #s *= remember_sign
    s *= -1
    fl.next('signal, phased')
    fl.plot(s)
    fl.next('signal, phased - image')
    fl.image(s)
    enhancement = s['t2':(-1e3,1e3)].C
    #enhancement = s.C
    enhancement.sum('t2').real
    enhanced = enhancement.data
    enhanced /= max(enhanced)
    fl.next(r'150$\mu$M TEMPOL E(p)')
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
    power_axis_W = r_[0,power_axis_W]
    fl.plot(power_axis_W[:-3],enhanced[:-3],'.',human_units=False)
    fl.plot(power_axis_W[-3:],enhanced[-3:],'o',human_units=False)
    xlabel('Power (W)')
    ylabel('Enhancement')
    if k_sigma:
        T1s = r_[1.897,1.915,2.551,3.217]
        p_dBm = r_[23.44,30.81,34.81]
        p_W = 10**(p_dBm/10.) * 1e-3
        p_W = r_[0,p_W]
        f = interpolate.interp1d(p_W,T1s)
        xnew = linspace(p_W[0],p_W[-1],25)
        ynew = f(xnew)
        figure()
        T1_p = ndshape([len(p_W)],['powers']).alloc(complex128)
        T1_p.setaxis('powers',p_W)
        T1_p['powers',:] = T1s[:]
        fl.next('T1 plot')
        fl.plot(T1_p,'o')
        fl.plot(xnew,ynew,'.')
        ks_p = ndshape([len(power_axis_W[:-3])],['powers']).alloc()
        ks_p.setaxis('powers',power_axis_W[:-3])
        ks_p['powers',:] = 1-enhanced[:-3]
        ks_p /= (395e-6*(9.822555e9/14.898292e6)*ynew)
        fl.next(r'395$\mu$M TEMPOL $k_{\sigma}$s(p)')
        fl.plot(ks_p,'.')
fl.show();quit()
