from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
fl = figlist_var()
for date,id_string in [
        ('191122','CPMG_2')
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
    nScans = s.get_prop('acq_params')['nScans']
    s.set_units('t','s')
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
    s.setaxis('nScans',r_[0:nScans])
    print ndshape(s)
    s.reorder('nScans',first=True)
    print ndshape(s)
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
    #s.mean('nScans',return_error=False)
    echo_center = ones_like(s.getaxis('nScans'),dtype=float)
    for x in xrange(nScans):
        echo_center[x] = abs(s['nScans',x])['tE',0].argmax('t2').data.item()
    s.setaxis('t2', lambda x: x-echo_center.mean())
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    fl.next('check center')
    fl.image(s)
    s.ft('t2')
    fl.next('before phased - real ft')
    fl.image(s.real)
    fl.next('before phased - imag ft')
    fl.image(s.imag)
    f_axis = s.fromaxis('t2')
    def costfun(p):
        zeroorder_rad,firstorder = p
        phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
        phshift *= exp(-1j*2*pi*zeroorder_rad)
        corr_test = phshift * s
        return (abs(corr_test.data.imag)**2)[:].sum()
    iteration = 0
    def print_fun(x, f, accepted):
        global iteration
        iteration += 1
        print (iteration, x, f, int(accepted))
        return
    sol = basinhopping(costfun, r_[0.,0.],
            minimizer_kwargs={"method":'L-BFGS-B'},
            callback=print_fun,
            stepsize=100.,
            niter=10,
            T=1000.
            )
    zeroorder_rad, firstorder = sol.x
    phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(-1j*2*pi*zeroorder_rad)
    s *= phshift
    print "RELATIVE PHASE SHIFT WAS {:0.1f}\us and {:0.1f}$^\circ$".format(
            firstorder,angle(zeroorder_rad)/pi*180)
    if s['nEchoes',0].data[:].sum().real < 0:
        s *= -1
    print ndshape(s)
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)
    s.ift('t2')
    fl.next('after phased - real')
    fl.image(s.real)
    fl.next('after phased - imag')
    fl.image(s.imag)
    s.rename('nEchoes','tE').setaxis('tE',tE_axis)
    data = s.C.sum('t2')
    fl.next('Fit decay')
    x = tE_axis
    print x
    print shape(x)
    print type(x)
    for scan_num in xrange(nScans):
        ydata = data['nScans',scan_num].data.real
        ydata /= max(ydata)
        print type(ydata)
        print shape(ydata)
        fl.plot(x,ydata, '.', alpha=0.4, label='data', human_units=False)
        fitfunc = lambda p, x: exp(-x/p[0])
        errfunc = lambda p_arg, x_arg, y_arg: fitfunc(p_arg, x_arg) - y_arg
        p0 = [0.2]
        p1, success = leastsq(errfunc, p0[:], args=(x, ydata))
        x_fit = linspace(x.min(),x.max(),5000)
        fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T2 = %0.2f ms)'%(p1[0]*1e3), human_units=False)
        xlabel('t (sec)')
        ylabel('Intensity')
        T2 = p1[0]
        print "T2:",T2,"s"
    fl.show();quit()
    s.ift('t2')
    even_echo_center = abs(s)['ph1',1]['tE',0].argmax('t2').data.item()
    odd_echo_center = abs(s)['ph1',-1]['tE',1].argmax('t2').data.item()
    print "EVEN ECHO CENTER:",even_echo_center,"s"
    print "ODD ECHO CENTER:",odd_echo_center,"s"
    s.setaxis('t2',lambda x: x-even_echo_center)
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    fl.next('check center before interleaving')
    fl.image(s)
    interleaved = ndshape(s)
    interleaved['ph1'] = 2
    interleaved['nEchoes'] /= 2
    interleaved = interleaved.rename('ph1','evenodd').alloc()
    #interleaved.copy_props(s).setaxis('t2',s.getaxis('t2').copy()).set_units('t2',s.get_units('t2'))
    interleaved.setaxis('t2',s.getaxis('t2').copy()).set_units('t2',s.get_units('t2'))
    interleaved.ft('t2',shift=True)
    interleaved.ift('t2')
    interleaved['evenodd',0] = s['ph1',1]['nEchoes',0::2].C.run(conj)['t2',::-1]
    interleaved['evenodd',1] = s['ph1',-1]['nEchoes',1::2]
    interleaved.ft('t2')
    fl.next('even and odd')
    fl.image(interleaved)
    phdiff = interleaved['evenodd',1]/interleaved['evenodd',0]*abs(interleaved['evenodd',0])
    fl.next('phdiff')
    fl.image(phdiff)
    phdiff *= abs(interleaved['evenodd',1])
    f_axis = interleaved.fromaxis('t2')
    def costfun(firstorder):
        phshift = exp(-1j*2*pi*f_axis*firstorder)
        return -1*abs((phdiff * phshift).data[:].sum())
    sol = minimize(costfun, ([0],),
            method='L-BFGS-B',
            bounds=((-1e-3,1e-3),)
            )
    firstorder = sol.x[0]
    phshift = exp(-1j*2*pi*f_axis*firstorder)
    phdiff_corr = phdiff.C
    phdiff_corr *= phshift
    zeroorder = phdiff_corr.data[:].sum().conj()
    zeroorder /= abs(zeroorder)
    fl.next('phdiff -- corrected')
    fl.image(phdiff_corr)
    print "Relative phase shift (for interleaving) was "        "{:0.1f}\us and {:0.1f}$^\circ$".format(
                firstorder/1e-6,angle(zeroorder)/pi*180)
    interleaved['evenodd',1] *= zeroorder*phshift
    interleaved.smoosh(['nEchoes','evenodd'],noaxis=True).reorder('t2',first=False)
    interleaved.setaxis('nEchoes',r_[1:nEchoes+1])
    f_axis = interleaved.fromaxis('t2')
    def costfun(p):
        zeroorder_rad,firstorder = p
        phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
        phshift *= exp(-1j*2*pi*zeroorder_rad)
        corr_test = phshift * interleaved
        return (abs(corr_test.data.imag)**2)[:].sum()
    iteration = 0
    def print_fun(x, f, accepted):
        global iteration
        iteration += 1
        print (iteration, x, f, int(accepted))
        return
    sol = basinhopping(costfun, r_[0.,0.],
            minimizer_kwargs={"method":'L-BFGS-B'},
            callback=print_fun,
            stepsize=100.,
            niter=100,
            T=1000.
            )
    zeroorder_rad, firstorder = sol.x
    phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(-1j*2*pi*zeroorder_rad)
    interleaved *= phshift
    print "RELATIVE PHASE SHIFT WAS {:0.1f}\us and {:0.1f}$^\circ$".format(
            firstorder,angle(zeroorder_rad)/pi*180)
    if interleaved['nEchoes',0].data[:].sum().real < 0:
        interleaved *= -1
    print ndshape(interleaved)
    interleaved.reorder('t2',first=True)
    fl.next('phased echoes, real - ft')
    fl.plot(interleaved.real)
    fl.next('phased echoes, imag - ft')
    fl.plot(interleaved.imag)
    interleaved = interleaved['t2':(-4e3,4e3)].C
    fl.next('phased echoes, real')
    fl.plot(interleaved.real)
    fl.next('phased echoes, imag')
    fl.plot(interleaved.imag)
    interleaved.rename('nEchoes','tE').setaxis('tE',tE_axis)
    data = interleaved.C.sum('t2')
    fl.next('Fit decay')
    x = tE_axis 
    ydata = data.data.real
    ydata /= max(ydata)
    fl.plot(x,ydata, '.', alpha=0.4, label='data', human_units=False)
    fitfunc = lambda p, x: exp(-x/p[0])
    errfunc = lambda p_arg, x_arg, y_arg: fitfunc(p_arg, x_arg) - y_arg
    p0 = [0.2]
    p1, success = leastsq(errfunc, p0[:], args=(x, ydata))
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T2 = %0.2f ms)'%(p1[0]*1e3), human_units=False)
    xlabel('t (sec)')
    ylabel('Intensity')
    T2 = p1[0]
    print "T2:",T2,"s"
