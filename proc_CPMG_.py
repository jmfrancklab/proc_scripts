from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
fl = figlist_var()
for date,id_string in [
        ('191206','CPMG_1_0'),
        ('191206','CPMG_1_1'),
        ('191206','CPMG_1_2'),
        ('191206','CPMG_1_3'),
        ('191206','CPMG_1_4'),
        ('191206','CPMG_1_5'),
        ('191206','CPMG_1_6'),
        ('191206','CPMG_1_7'),
        ('191206','CPMG_1_8'),
        ('191206','CPMG_1_9'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    s.set_units('t','s')
    orig_t = s.getaxis('t')
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    transient_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank = s.get_prop('acq_params')['deblank_us']*1e-6
    acq_time_s = orig_t[nPoints]
    tau_s = s.get_prop('acq_params')['tau_us']*1e-6
    pad_s = s.get_prop('acq_params')['pad_us']*1e-6
    tE_s = 2.0*p90_s + transient_s + acq_time_s + pad_s
    print("ACQUISITION TIME:",acq_time_s,"s")
    print("TAU DELAY:",tau_s,"s")
    print("TWICE TAU:",2.0*tau_s,"s")
    print("ECHO TIME:",tE_s,"s")
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*tE_s
    s.setaxis('t',None)
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis)
    s.ft('t2', shift=True)
    s.ft(['ph1'])
    s.ift('t2')
    s = s['ph1',1].C
    s.reorder('t2',first=True)
    s.sum('t2')
    #s.smoosh(['tE','t2'])
    #s.setaxis('tE',r_[0:nPoints*nEchoes])
    fl.next('Coherence pathway: smooshed')
    fl.plot(abs(s),human_units=False,alpha=0.5,
            label='%s'%id_string)
fl.show();quit()
