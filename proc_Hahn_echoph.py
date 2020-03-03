from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
<<<<<<< HEAD
for date,id_string,label_str in [
        ('200302','alex_probe_w33_fullMW','waterloading 33'),
        #('191206','echo_TEMPOL_1','TEMPOL'),
=======
t2 = symbols('t2')
filter_bandwidth = 5e3
color_choice = True
for date,id_string,label_str,color_str in [
        ('200302','alex_probe_w33_noMW','microwaves off','blue'),
        ('200302','alex_probe_w33_fullMW','microwaves on','red'),
>>>>>>> 4b319462e8a789a7533a72c91c324a924241b6e3
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    print ndshape(s)
    s.reorder('t',first=True)
    t2_axis = s.getaxis('t')[0:2048]
    s.setaxis('t',None)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.setaxis('t2',t2_axis)
    s.setaxis('nScans',r_[0:nScans])
    s.reorder('t2',first=False)
<<<<<<< HEAD
    s.ft('t2',shift=True)
    fl.next('raw data, chunked')
    fl.image(abs(s))
    s.ft(['ph1','ph2'])
    fl.next('coherence')
    fl.image(abs(s))
    s = s['ph1',1]['ph2',0].C
    s.mean('nScans',return_error=False)
    fl.next('signal')
    fl.plot(abs(s),label=label_str)
    slice_f = (-5e3,5e3)
    s = s['t2':slice_f].C
    s.ift('t2')
    max_data = abs(s.data).max()
    print max_data
    pairs = s.contiguous(lambda x: abs(x) > max_data*0.5)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    s.setaxis('t2',lambda x: x-peak_location.mean())
    s.register_axis({'t2':0})
    max_shift = diff(peak_location).item()/2
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    s_ft = s_sliced.C
    fl.next('sliced')
    fl.plot(s_ft)
    shift_t = nddata(r_[-0.1:0.1:200j]*max_shift, 'shift')
    s_foropt = s.C
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
    print s_foropt.getaxis('t2')
    print s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]]
    if ndshape(s_foropt)['t2'] % 2 == 0:
        s_foropt = s_foropt['t2',:-1]
    assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2+1] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
    ph0 = s_foropt['t2':0.0]
    ph0 /= abs(ph0)
    s_foropt /= ph0
    s_foropt /= max(abs(s_foropt.getaxis('t2')))
    # }}}
    residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
    residual.reorder('shift')
    print ndshape(residual)
    fl.next('residual')
=======
    fl.next('raw data -- coherence channels')
    s.ft(['ph2','ph1'])
    fl.image(s)
    fl.next('filtered + rough centered data %s'%label_str)
    s.ft('t2', shift=True)
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.image(s)
    s.ft('t2')
    k = s.C
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1],shift_val=1)
    fl.next('hermitian test')
>>>>>>> 4b319462e8a789a7533a72c91c324a924241b6e3
    fl.plot(residual)
    minpoint = residual.argmin()
    best_shift = minpoint['shift']
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
<<<<<<< HEAD
    ph0 = s['t2':0.0]
    ph0 /= abs(ph0)
    s /= ph0
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    fl.next('time domain')
    fl.plot(s_sliced)
    s_sliced.ft('t2')
    fl.next('Spectrum FT')
    fl.plot(s_sliced.real, alpha=0.5, label='real - %s'%label_str)
    fl.plot(s_sliced.imag, alpha=0.5, label='imag - %s'%label_str)
=======
    s *= exp(-s.getaxis('t2')/40e-3)
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
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    s = s['t2':(0,None)]
    fl.next('phased - time')
    fl.plot(s['ph2',-2]['ph1',1])
    s.ft('t2')
    #s.convolve('t2',7)
    if label_str == 'microwaves on':
        s *= -1
    fl.next('')
    s.name('')
    if color_choice:
        fl.plot(s['ph2',-2]['ph1',1],label='%s'%label_str,c=color_str)
    else:
        fl.plot(s['ph2',-2]['ph1',1],label='%s'%label_str)
>>>>>>> 4b319462e8a789a7533a72c91c324a924241b6e3
fl.show();quit()
