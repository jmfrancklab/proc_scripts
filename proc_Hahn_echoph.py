# this version of the code taken from checking out the only commit on dec 12 2019
# to process, use 'conda activate py2' assuming you have a proper python2 environment setup
# then 'git checkout py2' in pyspecdata

from pyspecdata import * 
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
for date,id_string,label_str in [
        ('191212','echo_5','water-Ni'),
        #('191206','echo_TEMPOL_1','TEMPOL'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    print(s.get_prop('acq_params'))
    s.reorder('t',first=True)
    nPhaseSteps = 8.
    stop_pt = nPoints/nPhaseSteps
    t2_axis = s.getaxis('t')
    t2_axis = t2_axis[:int(stop_pt)]
    s.setaxis('t',None)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.setaxis('t2',t2_axis)
    s.setaxis('nScans',r_[0:nScans])
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    fl.next('raw data, chunked')
    fl.image(abs(s))
    s.ft(['ph1','ph2'])
    fl.next('coherence')
    fl.image(abs(s))
    s = s['ph1',1]['ph2',0].C
    s.mean('nScans')
    fl.next('signal')
    fl.plot(abs(s))
    slice_f = (-5e3,5e3)
    s = s['t2':slice_f].C
    s.ift('t2')
    max_data = abs(s.data).max()
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
    fl.next('residual')
    fl.plot(residual)
    minpoint = residual.argmin()
    best_shift = minpoint['shift']
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
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
fl.show();quit()
