from pyspecdata import *
from scipy.optimize import minimize
date = '200110'
nodename = 'nutation'
freq_window = (-100,200)
with figlist_var(filename='CPMG_data.pdf') as fl:
    for id_string in [
        'CPMG_calib_1',
        ]:
        fl.basename = id_string
        filename = date+'_'+id_string+'.h5'
        s = nddata_hdf5(filename+'/'+nodename,
                directory = getDATADIR(exp_type='test_equip'))
        # {{{ pull a bunch of parameters from data file
        SW_kHz = s.get_prop('acq_params')['SW_kHz']
        nPoints = s.get_prop('acq_params')['nPoints']
        nEchoes = s.get_prop('acq_params')['nEchoes']
        nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
        nScans = s.get_prop('acq_params')['nScans']
        p90_s = s.get_prop('acq_params')['p90_us']
        transient_s = s.get_prop('acq_params')['deadtime_us']*1e-6
        deblank = s.get_prop('acq_params')['deblank_us']*1e-6
        tau_s = s.get_prop('acq_params')['tau_us']*1e-6
        # }}}
        s.set_units('p_90','s')
        s.setaxis('nScans',r_[0:nScans])
        s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
        s.setaxis('ph1',r_[0.,2.]/4)
        s.ft('t2',shift=True)
        fl.next('image, raw')
        fl.image(s, interpolation='bilinear')
        s.ft(['ph1'])
        fl.next('image, all coherence channels')
        s.mean('nScans',return_error=False)
        s.reorder('ph1')
        fl.image(s, interpolation='bilinear')
        fl.next('image')
        s = s['ph1',1]['t2':(-500,500)]
        fl.image(s, interpolation='bilinear')
        fl.next('image abs')
        fl.image(abs(s), interpolation='bilinear')
        fl.next('time domain')
        s.ift('t2')
        print "start of time axis",s.getaxis('t2')[0]
        assert s.getaxis('t2')[0] == 0.0, "why doesn't your time axis start at zero??"
        t_last = s.getaxis('t2')[-1]
        s.setaxis('t2', lambda x: x-t_last/2)# assume echo is symmetric
        fl.image(s, interpolation='bilinear')# remember that we are
        #                           interpolating between relatively few pixes
        #                           on the x axis
        fl.next('check that echo seems centered on $t=0$')
        forplot = abs(s).sum('p_90').sum('tE')
        fl.plot(forplot, label='before correction')
        x_center = forplot.argmax('t2').item()
        s.setaxis('t2', lambda x: x-x_center)
        fl.plot(abs(s).sum('p_90').sum('tE'), label='after ad-hoc correction')
        fl.next('freq domain after centering echo')
        s.ft('t2')
        fl.image(s, interpolation='bilinear')
        fl.next('freq domain after centering echo')
        fl.image(abs(s), interpolation='bilinear')
        axvline(x=freq_window[0], color='r')
        axvline(x=freq_window[1], color='r')
        fl.next('sum along frequency dimension')
        s = s['t2':freq_window]
        s.sum('t2')
        s.setaxis('tE', lambda x: x+2*x_center)# first point is not at zero,
        #                        but halfway through the first echo *plus* one
        #                        tau period
        s.reorder('tE',first=False)
        fl.image(s)
        fl.next('phase spectra')
        phcorr = nddata(r_[0:pi:50j],'phcorr')
        signal = abs((s*exp(-1j*phcorr)).sum('tE').real).sum('p_90')
        signal_power = (abs(s)**2).sum('tE').sum('p_90')
        success = signal/signal_power
        fl.plot(success)
        phcorr = success.argmax()
        #fl.plot(phcorr,'o', human_units=False) # should be able to do this,
        #                                   but it refuses zero-d data --
        #                                   should just use the units and go
        #                                   ahead, but leave for figlist
        #                                   upgrade
        phcorr = phcorr.item()
        s *= exp(-1j*phcorr)
        if (s['p_90',:5].sum('tE').sum('p_90')).item().real < 0:
            s *= -1
        fl.next('signal phase points')
        fl.plot(s.data.ravel().real,s.data.ravel().imag,'k.',
                alpha=0.1,mec='none')
        gca().set_aspect('equal')
        fl.next('processed data as function of $t_{90}$: waterfall')
        s.reorder('tE') # the lines appear to be the first dimension for waterfall
        s.waterfall()
        fl.next('processed data as function of $t_{90}$')
        fl.plot(s)
        #fl.image(abs(s['t2':(0,100)]))
        ##s.sum('t2')
        #s.setaxis('tE',r_[0:nEchoes])
        #print p90_s
        #fl.next('Coherence pathway: smooshed')
        #abs(s['t2':(0,100)].sum('t2')).waterfall()
        #fl.show()
