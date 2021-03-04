from pyspecdata import *
fl = figlist_var()
for date,id_string in [
        ('200115','CPMG_DNP_1')
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    #{{{ pulling acq params
    s.setaxis('power',r_[
        0,dBm2power(array(s.get_prop('meter_powers'))+20)]
        ).set_units('power','W')
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
    #s.setaxis('nScans',r_[0:nScans])
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
    s = s['ph1',1].C
    fl.next(id_string+' image plot coherence-- ft ')
    fl.image(s)
    s.ift('t2')
    fl.next(id_string+' image plot coherence ')
    fl.image(s)
    #s.mean('nScans',return_error=False)
    s.reorder('t2',first=True)
    print(ndshape(s))
    t2_max = zeros_like(s.getaxis('power'))
    for x in range(len(s.getaxis('power'))):
        t2_max[x] = abs(s['power',x]['tE',0]).argmax('t2',raw_index=True).data
    s.setaxis('t2',lambda t: t -s.getaxis('t2')[int(t2_max.mean())])
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    print(ndshape(s))
    s.reorder('nEchoes',first=True)
    print(ndshape(s))
    s.ft('t2')
    # as of right now, this plot doesn't show anything meaningful
    fl.next('check center')
    fl.image(s)
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
    s *= phshift
    print("RELATIVE PHASE SHIFT WAS {:0.1f}\\us and {:0.1f}$^\circ$".format(
            firstorder,angle(zeroorder_rad)/pi*180))
    if s['nEchoes':0]['t2':0]['power',-4].item().real > 0:
        print(s['nEchoes':0]['t2':0]['power',-4].item().real)
        print("Sign correction")
        s *= -1
        print(s['nEchoes':0]['t2':0]['power',-4].item().real)
    print(ndshape(s))
    s.reorder('power',first=True)
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)
    T2_values = ones(ndshape(s)['power'])
    find_T2 = True
    #{{{ find T2
    if find_T2:
        for k in range(ndshape(s)['power']):
        #for k in [1,5,10,15,19]:
            data = s['t2':0]['power',k]
            fl.next('Echo decay (power = %f W)'%s.getaxis('power')[k])
            x = tE_axis
            ydata = data.data.real
            if ydata[0] < 1:
                ydata *= -1
            ydata /= max(ydata)
            fl.plot(x,ydata, '.', alpha=0.4, label='data', human_units=False)
            raise RuntimeEror("Need to set up fitdata using new capabilities")
            data.guess([0.1,100.0]) <-- should be a dictionary
            fl.plot(data.eval(400))
            T2_values[k] = T2
            print("T2:",T2,"s")
        fl.next('T2 vs power')
        fl.plot(s.getaxis('power')[:-3],T2_values[:-3],'.')
        xlabel('Power (W)')
        ylabel('T2 (seconds)')
    #}}}
    print(T2_values)
    fl.show();quit()
    enhancement = s['t2':0]['nEchoes',0].C
    enhanced = enhancement.data[1:].real
    enhanced /= max(enhanced)
    fl.next('150 uL TEMPOL enhancement (first echo of CPMG train)')
    fl.plot(s.getaxis('power')[:-3],enhanced[:-3],'.',human_units=False)
    #fl.plot(power_axis_W[-3:],enhanced[-3:],'o',human_units=False)
    xlabel('Power (W)')
    ylabel('Enhancement')

fl.show();quit()


