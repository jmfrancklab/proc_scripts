from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
fl = figlist_var()
for date,id_string,label_str in [
        ('200115','CPMG_26_3p0','3.0'),
        ('200115','CPMG_26_3p1','3.1'),
        ('200115','CPMG_26_3p2','3.2'),
        ('200115','CPMG_26_3p3','3.3'),
        ('200115','CPMG_26_3p4','3.4'),
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
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    s.set_units('t','s')
    twice_tau = deblank_s + 2*p90_s + deadtime_s + pad_start_s + acq_time_s + pad_end_s + marker_s
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*twice_tau
    s.setaxis('t',None)
    s.setaxis('nScans',r_[0:nScans])
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis)
    #fl.next(id_string+'raw data - chunking')
    #fl.image(s)
    s.ft('t2', shift=True)
    #fl.next(id_string+'raw data - chunking ft')
    #fl.image(s)
    s.ft(['ph1'])
    #fl.next(id_string+' image plot coherence-- ft ')
    #fl.image(s)
    s.ift('t2')
    s.reorder('nScans',first=True)
    #fl.next(id_string+' image plot coherence ')
    #fl.image(s, interpolation='bilinear')
    s = s['ph1',1].C
    s.mean('nScans',return_error=False)
    s.reorder('t2',first=True)
    echo_center = abs(s)['tE',0].argmax('t2').data.item()
    s.setaxis('t2', lambda x: x-echo_center)
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    fl.next('check center')
    fl.image(s)
    s.ft('t2')
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
            niter=100,
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
    #data = s['t2':0]
    data = s['t2':(-100,100)].sum('t2')
    fl.next('Echo decay')
    x = tE_axis
    ydata = data.data.real
    ydata /= max(ydata)
    fl.plot(x,ydata,'-o', alpha=0.7, label='%s'%label_str, human_units=False)
    fitfunc = lambda p, x: p[0]*exp(-x*p[1])
    errfunc = lambda p_arg, x_arg, y_arg: fitfunc(p_arg, x_arg) - y_arg
    p0 = [0.1,100.0,-3.0]
    p1, success = leastsq(errfunc, p0[:], args=(x, ydata))
    assert success == 1, "Fit did not succeed"
    T2 = 1./p1[1]
    print T2
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T2 = %0.2f ms)'%(T2*1e3), human_units=False)
    xlabel('t (sec)')
    ylabel('Intensity')
    print "T2:",T2,"s"
save_fig = False
if save_fig:
    savefig('20200108_CPMG_trials.png',
            transparent=True,
            bbox_inches='tight',
            pad_inches=0,
            legend=True)
fl.show()
