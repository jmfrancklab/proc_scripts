from pyspecdata import *
from scipy.optimize import minimize
fl = figlist_var()
date = '191206'
for id_string in [
    'CPMG_TEMPOL_calib_3',
    ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'nutation'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip' ))
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    orig_t = s.getaxis('t')
    p90_s = s.get_prop('acq_params')['p90_us']
    transient_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank = s.get_prop('acq_params')['deblank_us']*1e-6
    acq_time_s = orig_t[nPoints]
    tau_s = s.get_prop('acq_params')['tau_us']*1e-6
    pad_range = zeros_like(p90_s)
    for x in xrange(len(pad_range)):
        pad_range[x] = 2.0*tau_s - transient_s - acq_time_s - 2.0*p90_s[x]*1e-6 - deblank
    s.set_units('p_90','s')
    acq_time_s = orig_t[nPoints]
    t2_axis = linspace(0,acq_time_s,nPoints)
    s.setaxis('t',None)
    s.reorder('t',first=True)
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('t2',t2_axis)
    s.reorder('t2',first=False)
    s.reorder('p_90',first=True)
    s.ft('t2',shift=True)
    fl.next('image, raw')
    fl.image(s)
    s.ft(['ph1'])
    fl.next('image, all coherence channels')
    fl.image(s)
    s = s['ph1',1].C
    fl.next(id_string+'image')
    fl.image(s)
    fl.next(id_string+'image abs')
    fl.image(abs(s))
    s = s['t2':(-500,500)]
    s.sum('t2')
    s.reorder('tE',first=True)
    s.setaxis('tE',r_[0:nEchoes])
    fl.next('Coherence pathway: smooshed')
    abs(s).waterfall()
    fl.show()
fl.show();quit()
