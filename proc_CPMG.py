from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
fl = figlist_var()
for date,id_string in [
        ('191206','CPMG_2')
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
    print ndshape(s)
    fl.next(id_string+'raw data ')
    fl.plot(s.real,alpha=0.4)
    fl.plot(s.imag,alpha=0.4)
    fl.plot(abs(s),':',c='k',alpha=0.4)
    orig_t = s.getaxis('t')
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    transient_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank = s.get_prop('acq_params')['deblank_us']*1e-6
    acq_time_s = orig_t[nPoints]
    tau_s = s.get_prop('acq_params')['tau_us']*1e-6
    pad_s = s.get_prop('acq_params')['pad_us']*1e-6
    tE_s = 2.0*p90_s + transient_s + acq_time_s + pad_s
    print "ACQUISITION TIME:",acq_time_s,"s"
    print "TAU DELAY:",tau_s,"s"
    print "TWICE TAU:",2.0*tau_s,"s"
    print "ECHO TIME:",tE_s,"s"
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*tE_s
    s.setaxis('t',None)
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis)
    fl.next(id_string+'raw data - chunking')
    fl.image(s)
    s.ft('t2', shift=True)
    fl.next(id_string+'raw data - chunking ft')
    fl.image(s)
    s.ft(['ph1'])
    fl.next(id_string+' image plot coherence-- ft ')
    fl.image(s)
    s.ift('t2')
    fl.next(id_string+' image plot coherence ')
    fl.image(s)
    s = s['ph1',1].C
    s.reorder('t2',first=True)
    s.ft('t2')
    slice_f = (-4e3,4e3)
    s = s['t2':slice_f].C
    s.ift('t2')
    first_s = s['tE',0].C
    max_data = abs(first_s.data).max()
    pairs = first_s.contiguous(lambda x: abs(x) > max_data*0.5)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    print peak_location
    s.setaxis('t2', lambda x: x-peak_location.mean())
    first_s.setaxis('t2', lambda x: x-peak_location.mean())
    s.register_axis({'t2':0})
    first_s.register_axis({'t2':0})
    max_shift = diff(peak_location).item()/2
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    s_ft = s_sliced.C
    shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
    t2_decay = exp(-first_s.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
    s_foropt = first_s.C
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    s_foropt /= t2_decay
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
    residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
    residual.reorder('shift')
    print ndshape(residual)
    print residual
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
    #fl.next('time domain')
    #fl.plot(s_sliced)
    s.ft('t2')
    s_sliced.ft('t2')
    fl.next('Spectrum FT')
    fl.plot(s_sliced.real, alpha=0.5)
    fl.next('Spectrum FT - imag')
    fl.plot(s_sliced.imag, alpha=0.5)
    data = s_sliced['t2':(-0.25e3,0.4e3)].C
    ydata = data.C.real.sum('t2')
    ydata = ydata.data
    fl.next('Fit decay')
    x = tE_axis 
    ydata /= max(ydata)
    fl.plot(x,ydata, '.', alpha=0.4, label='data', human_units=False)
    #fl.next('after phased - real ft')
    #fl.image(s_sliced.real)
    #fl.next('after phased')
    #s_sliced.reorder('t2',first=True)
    #fl.plot(s_sliced.real)
    #fl.next('after phased - imag ft')
    #fl.image(s_sliced.imag)
    #s_sliced.ift('t2')
    #fl.next('after phased - real')
    #fl.image(s_sliced.real)
    #fl.next('after phased - imag')
    #fl.image(s_sliced.imag)
    s_sliced.setaxis('tE',tE_axis)
    #data = s_sliced['t2':(-0.1e3,0.2e3)].C.sum('t2')
    fitfunc = lambda p, x: exp(-x/p[0])
    errfunc = lambda p_arg, x_arg, y_arg: fitfunc(p_arg, x_arg) - y_arg
    p0 = [1.0]
    p1, success = leastsq(errfunc, p0[:], args=(x, ydata))
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T2 = %0.2f ms)'%(p1[0]*1e3), human_units=False)
    xlabel('t (sec)')
    ylabel('Intensity')
    T2 = p1[0]
    print "T2:",T2,"s"
fl.show()
