from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
def get_W(dBm):
    return 10**(dBm/10.)*1e-3
dBm_list = [0., 30., 34., 36.]
W_list = ones_like(dBm_list)
for x in range(4):
    W_list[x] = get_W(dBm_list[x])
enhancement = []
find_phase_params = True # phase params found for first dataset will be applied
                         # to all subsequently processed datasets
for date,id_string,label_string in [
        ('191031','echo_5_4','no microwaves'),
        ('191031','echo_5_mw_30dBm','+30 dBm microwaves'),
        ('191031','echo_5_mw_34dBm','+34 dBm microwaves'),
        ('191031','echo_5_mw_36dBm_2','+36 dBm microwaves'),
        ]:
    title_string = 'unenhanced'
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    print(ndshape(s))
    print(s.get_prop('acq_params'))
    quit()
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
    s = s['ph1',1]['ph2',0]
    s.ft('t2', shift=True)
    fl.next('frequency domain')
    fl.plot(abs(s), human_units=False)# no human units to make sure
    slice_f = (-1e3,1e3)
    s = s['t2':slice_f]
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
    if find_phase_params:
        shift_t = nddata(r_[-1:1:1000j]*max_shift, 'shift')
        t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:1000j],'R2'))
        s_foropt = s.C
        s_foropt.ft('t2')
        s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
        s_foropt.ift('t2')
        s_foropt /= t2_decay
        s_foropt = s_foropt['t2':(-max_shift,max_shift)]
        print(s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]])
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
        find_phase_params = False
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    ph0 = s['t2':0.0]
    ph0 /= abs(ph0)
    s /= ph0
    fl.next('Aer spectra')
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    fl.plot(s_sliced.real, alpha=0.5, label='%s'%label_string)
    enhancement.append(s_sliced.real.sum('t2').item())
fl.next('E(p)')
fl.plot(W_list,array(enhancement),'o-')
fl.show()
    
