from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
for date,id_string,label_string in [
        ('191031','echo_5_4','no microwaves'),
        ]:
    title_string = 'unenhanced'
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    print ndshape(s)
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1]).set_units('t2','s')
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.reorder('t2',first=False)
    s.ft(['ph1','ph2'])
    fl.next(id_string+'raw data - chunking coh')
    fl.image(s)
    s = s['ph1',1]['ph2',0]
    fl.next(id_string+'select pathway')
    fl.plot(s)
    s.ft('t2', shift=True)
    fl.next('frequency domain')
    fl.plot(abs(s), human_units=False)# no human units to make sure
    slice_f = (-1e3,1e3)
    s = s['t2':slice_f]
    s.ift('t2')
    fl.next('frequency filtered')
    fl.plot(s,human_units=False)
    max_data = abs(s.data).max()
    pairs = s.contiguous(lambda x: abs(x) > max_data*0.5)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    s.setaxis('t2',lambda x: x-peak_location.mean())
    s.register_axis({'t2':0})
    fl.next('crude centering')
    fl.plot(s)
    max_shift = diff(peak_location).item()/2
    fl.next('spectrum with crude centering')
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    s_ft = s_sliced.C
    fl.plot(s_sliced)
    shift_t = nddata(r_[-1:1:1000j]*max_shift, 'shift')
    t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:1000j],'R2'))
    s_foropt = s.C
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    s_foropt /= t2_decay
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
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
    fl.next('cost function')
    residual.reorder('shift')
    fl.image(residual)
    fl.plot(residual.C.argmin('shift').name('shift'),'x')
    minpoint = residual.argmin()
    best_shift = minpoint['shift']
    best_R2 = minpoint['R2']
    fl.plot(best_R2,best_shift,'o')
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    ph0 = s['t2':0.0]
    ph0 /= abs(ph0)
    s /= ph0
    fl.next('spectrum with optimized centering')
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    fl.plot(s_ft.real, alpha=0.5, label='pre-herm, real')
    fl.plot(s_ft.imag, alpha=0.5, label='pre-herm,imag')
    fl.plot(s_sliced.real, alpha=0.5, label='post-herm, real')
    fl.plot(s_sliced.imag, alpha=0.5, label='post-herm, imag')
fl.show()
    
