from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
rcParams["savefig.transparent"] = True
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 5e3
for date,id_string in [
    #('200122','echo_DNP_TCM51C_3'),
    #('200127','echo_DNP_TCM51C_1'),
    #('200128','echo_DNP_TCM118C_1'),
    #('200130','echo_DNP_1'),
    #('200130','echo_DNP_2'),
    #('200130','echo_DNP_3'),
    #('191118','echo_DNP_3'),
    #('191217','echo_DNP_1'),
    #('200130','echo_DNP_5'),
    #('200130','echo_DNP_AG'),
    #('200225','DNP_echo_1'),
    #('200221','DNP_S179R1apR_one'),
    ('200302','DNP_echo_w33_1'),
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
    #s = s['ph1',1]['ph2',0].C
    fl.next('viz - signal')
    fl.image(s)
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.next('Centering')
    fl.image(s)
    residual,best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1], shift_val=6.0)
    fl.next('hermitian test')
    fl.plot(residual)
    print("best shift is",best_shift)
    # {{{ slice out the FID appropriately and phase correct
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    fl.next('time domain after hermitian test')
    fl.image(s)
    ph0 = s['t2':0]['ph2',-2]['ph1',1]
    print(ndshape(ph0))
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0, fl=fl)
        print('phasing dimension as one')
    else:
        print("there is only one dimension left -- standard 1D zeroth order phasing")
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('time domain -- after hermitian function test and phasing')
    fl.plot(s['ph2',-2]['ph1',1])
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    fl.next('FID slice, time domain -- after hermitian function test and phasing')
    fl.plot(s['ph2',-2]['ph1',1])
    fl.next('FID slice, freq domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.plot(s['ph2',-2]['ph1',1])
    fl.show();quit()
    s = s['ph2',-2]['ph1',1]
    s *= -1
    enhancement = s['t2':(-0.1e3,0.2e3)].C
    #enhancement = s.C
    enhancement.sum('t2').real
    fl.next('plot')
    fl.plot(enhancement,'.')
    enhanced = enhancement.data
    enhanced /= max(enhanced)
    #fl.next(r'TEMPOL - Enhancement vs. Power')
    fl.next(r'')
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
    power_axis_W = r_[0,power_axis_W]
    fl.plot(power_axis_W[:-3],enhanced[:-3],'o',c='k',human_units=False,label='%s'%date)
    fl.plot(power_axis_W[-3:],enhanced[-3:],'x',c='red',human_units=False,label='%s'%date)
    xlabel('Power \ W')
    ylabel('Enhancement')
fl.show();quit()
