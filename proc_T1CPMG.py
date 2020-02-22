from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]

for date,id_string in [
        ('200221','T1CPMG_TEMPOLgel_1')
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
            #{{{ pulling acq params
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    nScans = s.get_prop('acq_params')['nScans']
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    deadtime_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank_s = s.get_prop('acq_params')['deblank_us']*1e-6
    marker_s = s.get_prop('acq_params')['marker_us']*1e-6
    tau1_s = s.get_prop('acq_params')['tau1_us']*1e-6
    pad_start_s = s.get_prop('acq_params')['pad_start_us']*1e-6
    pad_end_s = s.get_prop('acq_params')['pad_end_us']*1e-6
    #}}}
    s.set_units('t','s')
    fl.next('raw data - no clock correction')
    fl.image(s)
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    s.set_units('t','s')
    twice_tau = deblank_s + 2*p90_s + deadtime_s + pad_start_s + acq_time_s + pad_end_s + marker_s
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*twice_tau
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis).set_units('t2','s')
    vd_list = s.getaxis('vd')
    s.setaxis('vd',vd_list).set_units('vd','s')
    s.ft('t2',shift=True)
    #{{{ for applying clock correction
    clock_corr = True
    if clock_corr:
        clock_correction = 1.785
        s *= exp(-1j*s.fromaxis('vd')*clock_correction)
        fl.next('raw data - clock correction')
        fl.image(s['tE',0])
    #}}}
    s.ift('t2')
    s.ft(['ph1'])
    s = s['ph1',1].C
    fl.next(id_string+' image plot coherence')
    fl.image(s)
    s.ft('t2')
    fl.next(id_string+' image plot coherence -- ft')
    fl.image(s)
    s.ift('t2')
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    echo_center = abs(s)['nEchoes',0]['vd',0].argmax('t2').data.item()
    s.setaxis('t2', lambda x: x-echo_center)
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
        corr_test = phshift * s['vd',-1]
        return (abs(corr_test.data.imag)**2)[:].sum()
    iteration = 0
    def print_fun(x, f, accepted):
        global iteration
        iteration += 1
        print((iteration, x, f, int(accepted)))
        return
    sol = basinhopping(costfun, r_[0.,0.],
            minimizer_kwargs={"method":'L-BFGS-B'},
            callback=print_fun,
            stepsize=100.,
            niter=100,
            T=1000.
            )
    zeroorder_rad, firstorder = sol.x
    phshift = exp(1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(1j*2*pi*zeroorder_rad)
    s *= phshift
    print("RELATIVE PHASE SHIFT WAS {:0.1f}\\us and {:0.1f}$^\circ$".format(
            firstorder,angle(zeroorder_rad)/pi*180))
    #if s['nEchoes',0].data[:].sum().real < 0:
    #    s *= -1
    print(ndshape(s))
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)
    fl.show();quit()
    s.ift('t2')
    fl.next('after phased - real')
    fl.image(s.real)
    fl.next('after phased - imag')
    fl.image(s.imag)
    fl.next('echoes')
    view_echoes = False
    if view_echoes:
        for x,y in enumerate(vd_list):
            fl.plot(s['nEchoes',0]['vd',x],label='%s'%x)
        fl.next('echoes imag')
        for x,y in enumerate(vd_list):
            fl.plot(s.imag['nEchoes',0]['vd',x],label='%s'%x)
    ##
    fl.show();quit()
    even_echo_center = abs(s)['ph1',1]['vd',0]['nEchoes',0].argmax('t2').data.item()
    odd_echo_center = abs(s)['ph1',-1]['vd',0]['nEchoes',1].argmax('t2').data.item()
    print("EVEN ECHO CENTER:",even_echo_center,"s")
    print("ODD ECHO CENTER: ",odd_echo_center,"s")
    s.setaxis('t2',lambda x: x-even_echo_center)
    fl.next('check center before interleaving')
    fl.image(s)
    interleaved = ndshape(s)
    interleaved['ph1'] = 2
    interleaved['nEchoes'] /= 2
    interleaved = interleaved.rename('ph1','evenodd').alloc()
    interleaved.copy_props(s).setaxis('t2',s.getaxis('t2').copy()).set_units('t2',s.get_units('t2'))
    interleaved['evenodd',0] = s['ph1',1]['nEchoes',0::2].C.run(conj)['t2',::-1]
    interleaved['evenodd',1] = s['ph1',-1]['nEchoes',1::2]
    interleaved.ft('t2')
    fl.next('even and odd')
    fl.image(interleaved)
    phdiff = interleaved['evenodd',1]/interleaved['evenodd',0]*abs(interleaved['evenodd',0])
    fl.next('phdiff')
    fl.image(phdiff)
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
    print("Relative phase shift (for interleaving) was "        "{:0.1f}\\us and {:0.1f}$^\circ$".format(
                firstorder/1e-6,angle(zeroorder)/pi*180))
    interleaved['evenodd',1] *= zeroorder*phshift
    interleaved.smoosh(['nEchoes','evenodd'],noaxis=True).reorder('t2',first=False)
    interleaved.setaxis('nEchoes',r_[1:nEchoes+1])
    interleaved.ift('t2')
    fl.next('interleaved')
    fl.image(interleaved)
    interleaved.ft('t2')
    fl.next('interleaved -- ft')
    fl.image(interleaved)
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
        print((iteration, x, f, int(accepted)))
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
    print("RELATIVE PHASE SHIFT WAS {:0.1f}\\us and {:0.1f}$^\circ$".format(
            firstorder,angle(zeroorder_rad)/pi*180))
    if interleaved['nEchoes',0]['vd',0].data[:].sum().real > 0:
        interleaved *= -1
    fl.next('interleaved -- phased ft')
    fl.image(interleaved)
    fl.next('final real data ft')
    fl.image(interleaved.real)
    fl.next('final imag data ft')
    fl.image(interleaved.imag)
    interleaved = interleaved['t2':(-4e3,4e3)].C
    interleaved.ift('t2')
    print(ndshape(interleaved))
    interleaved.rename('nEchoes','tE').setaxis('tE',tE_axis)
    #for index,val in enumerate(vd_list):
    for index in r_[0:len(vd_list):1]:
        data_T2 = interleaved['vd',index].C.sum('t2')
        x = tE_axis 
        ydata = data_T2.data.real
        if ydata.sum() < 0:
            ydata *= -1
        ydata /= max(ydata)
        fl.next('T2 decay along vd: data')
        fl.plot(x,ydata, '.', alpha=0.4, label='data %d'%index, human_units=False)
        fl.next('Fit decay')
        fl.plot(x,ydata, '.', alpha=0.4, label='data %d'%index, human_units=False)
        fitfunc = lambda p, x: exp(-x/p[0])
        errfunc = lambda p_arg, x_arg, y_arg: fitfunc(p_arg, x_arg) - y_arg
        p0 = [0.2]
        p1, success = leastsq(errfunc, p0[:], args=(x, ydata))
        x_fit = linspace(x.min(),x.max(),5000)
        fl.next('T2 decay along vd: fit')
        fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit %d (T2 = %0.2f ms)'%(index,p1[0]*1e3), human_units=False)
        fl.next('Fit decay')
        fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit %d (T2 = %0.2f ms)'%(index,p1[0]*1e3), human_units=False)
        xlabel('t (sec)')
        ylabel('Intensity')
        T2 = p1[0]
        print("T2:",T2,"s")
fl.show();quit()
