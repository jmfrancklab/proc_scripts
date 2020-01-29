from pyspecdata import *
from sympy import symbols
transients_not_averaged = True
calculate_enhancements = True
import re
t2 = symbols('t2')
slice_f = (-1e3,1e3)
exp_tuples = [
        ('191111','echo_4','gradient off'),
        ('191111','echo_4_2','gradient off'),
        ('191111','echo_4_3','gradient off'),
        ('191111','echo_4_on','gradient on'),
        ('191111','echo_4_on_2','gradient on'),
        ('191111','echo_4_on_3','gradient on'),
        ]
if calculate_enhancements:
    dBm_re = re.compile(r'\b([0-9]+) *dBm\b')
    dBm_gen = (float(dBm_re.search(j).groups()[0])
            if dBm_re.search(j) is not None else -999
            for _,_,j in exp_tuples)
    def get_W(dBm):
        return 10**(dBm/10.)*1e-3
    W_list = map(get_W,dBm_gen)
    enhancement = ndshape([len(exp_tuples)],['power']).alloc().setaxis('power',
            array(list(W_list)))
with figlist_var() as fl:
    for j,(date,id_string,label_str) in enumerate(exp_tuples):
        filename = date+'_'+id_string+'.h5'
        nodename = 'signal'
        s = nddata_hdf5(filename+'/'+nodename,
                directory = getDATADIR(
                    exp_type = 'test_equip'))
        nPoints = s.get_prop('acq_params')['nPoints']
        nScans = s.get_prop('acq_params')['nScans']
        SW_kHz = s.get_prop('acq_params')['SW_kHz']
        nEchoes = s.get_prop('acq_params')['nEchoes']
        nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1]).set_units('t2','s')
        s.setaxis('ph2',r_[0.,2.]/4)
        s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
        if nScans > 1 and transients_not_averaged:
            s.setaxis('nScans',r_[0:nScans])
        s_shape = ndshape(s)
        fl.text(r'this data has %d phase cycle steps, it has been \textbf{averaged %f times} its \textbf{dw=%f \us} and its \textbf{aq=%f s}'%(nPhaseSteps,nScans,1e3/SW_kHz,s_shape['t2']*1e-3/SW_kHz))
        assert s_shape['t2'] == nPoints, "data doesn't appear to be the right shape: "+str(s_shape)
        rough_center = abs(s).mean_all_but('t2').argmax('t2').item()
        s.setaxis(t2-rough_center)
        s.ft('t2',shift=True).reorder(['ph2','ph1'])
        fl.next('raw data -- FT')
        fl.image(s['t2':(-300,300)])
        fl.next('raw data -- time domain')
        s.ift('t2')
        fl.image(s)
        s.ft(['ph1','ph2'])
        s.ft('t2')
        fl.next('coherence levels')
        fl.image(s)
        s = s['ph1',1]['ph2',-2]
        s = s['t2':slice_f]
        if nScans > 1:
            s.mean('nScans',return_error=False)
        fl.next('overlayed spectra')
        fl.plot(s)
        s.ift('t2')
        max_data = abs(s.data).max()
        pairs = s.contiguous(lambda x: abs(x) > max_data*0.5)
        longest_pair = diff(pairs).argmax()
        peak_location = pairs[longest_pair,:]
        s.setaxis('t2',lambda x: x-peak_location.mean())
        s.register_axis({'t2':0})
        max_shift = diff(peak_location).item()/2
        shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
        t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
        s_foropt = s.C
        s_foropt.ft('t2')
        s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
        s_foropt.ift('t2')
        s_foropt /= t2_decay
        s_foropt = s_foropt['t2':(-max_shift,max_shift)]
        logger.debug(s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]])
        if ndshape(s_foropt)['t2'] % 2 == 0:
            s_foropt = s_foropt['t2',:-1]
        assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2+1] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
        ph0 = s_foropt['t2':0.0]
        ph0 /= abs(ph0)
        s_foropt /= ph0
        s_foropt /= max(abs(s_foropt.getaxis('t2'))) # this is not needed?
        # }}}
        residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
        residual.reorder('shift')
        minpoint = residual.argmin()
        best_shift = minpoint['shift']
        best_R2 = minpoint['R2']
        s.ft('t2')
        s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
        s.ift('t2')
        ph0 = s['t2':0.0]
        ph0 /= abs(ph0)
        s /= ph0
        fl.next('Spectra with Hermitian phasing')
        s_sliced = s['t2':(0,None)].C
        s_sliced['t2',0] *= 0.5
        s_sliced.ft('t2')
        fl.plot(s_sliced.real, alpha=0.5, label='%s'%label_str)
        if calculate_enhancements:
            # this can be improved -- see notes above
            enhancement['power',j] = s_sliced.real.sum('t2').item()
    if calculate_enhancements:
        fl.next('E(p)')
        fl.plot(enhancement,'o-')
