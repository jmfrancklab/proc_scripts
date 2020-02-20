from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
for date,id_string,label_str in [
        ('200219','AG_probe_1','n'),
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
    t2_axis = s.getaxis('t')[0:nPoints/nPhaseSteps]
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.setaxis('nScans',r_[0:nScans])
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    fl.next('raw data, chunked')
    fl.image(abs(s)['t2':(-250,250)])
    s.ft(['ph1','ph2'])
    fl.next('coherence')
    fl.image(abs(s)['t2':(-250,250)])
    s = s['ph1',1]['ph2',-2].C
    s.mean('nScans')
    #s.mean('nScans',return_error=False)
    fl.next('plotting selected coherence channel')
    fl.plot(s.real, alpha=0.4, label='%s'%label_str)
    slice_f = (-3e3,3e3)
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
    s_ft = s_sliced.C
    fl.next('sliced')
    fl.plot(s_ft)
    shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
    t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
    s_foropt = s.C
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    s_foropt /= t2_decay
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
    print s_foropt.getaxis('t2')
    #print s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]]
    if ndshape(s_foropt)['t2'] % 2 == 0:
        s_foropt = s_foropt['t2',:-1]
    #assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2+1] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
    ph0 = s_foropt['t2':0.0]
    ph0 /= abs(ph0)
    s_foropt /= ph0
    s_foropt /= max(abs(s_foropt.getaxis('t2')))
    # }}}
    residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
    residual.reorder('shift')
    print ndshape(residual)
    minpoint = residual.argmin()
    best_shift = minpoint['shift']
    best_R2 = minpoint['R2']
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    ph0 = s['t2':0.0]
    ph0 /= abs(ph0)
    s /= ph0
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    fl.next('Spectrum - freq domain')
    fl.plot(s_sliced.real, alpha=0.5, label='%s'%filename)
fl.show();quit()
