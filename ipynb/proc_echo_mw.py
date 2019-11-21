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
    s.reorder('t2',first=True)
    t2_max = zeros_like(s.getaxis('power'))
    for x in xrange(len(s.getaxis('power'))):
        t2_max[x] = abs(s['power',x]['ph1',1]['ph2',0]).argmax('t2',raw_index=True).data
    print t2_max
    s.setaxis('t2',lambda t: t - s.getaxis('t2')[int(t2_max.mean())])
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2',shift=True)
    s.reorder('t2',first=False)
    fl.next('FID at t=0, then FT')
    fl.image(s)
    s.ft(['ph1','ph2'])
    fl.next('Coherence levels')
    fl.image(s)
fl.show()
