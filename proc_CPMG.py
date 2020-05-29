from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import postproc_dict 
fl = figlist_var()

for searchstr, exp_type, nodename, postproc, label_str in [
        #('200221_CPMG_TEMPOLgel_3p0_1','deadtime=5'),
        #('200221_CPMG_TEMPOLgel_2p9_1','deadtime=5'),
        ('200304_CPMG_2p6_1','test_equip','signal','CPMG','deadtime=5'),
        #('200305_CPMG_3p5_2','deadtime=5'),
        #('200305_CPMG_3p6_2','deadtime=5'),
        #('200305_CPMG_3p7_2','test_equip','signal','CPMG','deadtime=5'),
        #('200305_CPMG_3p7_3','deadtime=5'),
        #('200305_CPMG_3p8_2','deadtime=5'),
        #('200305_CPMG_3p9_2','deadtime=5'),
        #('200305_CPMG_4p0_1','deadtime=5'),
        ]:
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename, postproc=postproc, lookup=postproc_dict)
    logging.basicConfig(filename="logfile.txt")
    stderrLogger=logging.StreamHandler()
    stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logging.getLogger().addHandler(stderrLogger)
    nEchoes = s.get_prop('acq_params')['nEchoes']
    fl.next('raw data - chunking ft')
    fl.image(s)
    s.ft(['ph1'])
    fl.next(' image plot coherence-- ft ')
    fl.image(s)
    s.ift('t2')
    s.reorder('nScans',first=True)
    s = s['ph1',1].C
    s.mean('nScans')
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
        logger.info(strm((iteration, x, f, int(accepted))))
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
    logging = init_logging("info")
    logging.info("RELATIVE PHASE SHIFT WAS %0.1f\\us and %0.1f$^\circ$", firstorder, angle(zeroorder_rad)/pi*180)
    quit()
    if s['nEchoes',0].data[:].sum().real < 0:
        s *= -1
    logger.info(strm(ndshape(s)))
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)
    data = s['t2':(-200,200)].sum('t2')
    fl.next('Echo decay')
    x = s.getaxis('nEchoes')
    ydata = data.data.real
    ydata /= max(ydata)
    fl.plot(x,ydata,'-o', alpha=0.7, label='%s'%label_str, human_units=False)
    fitfunc = lambda p, x: p[0]*exp(-x*p[1])
    errfunc = lambda p_arg, x_arg, y_arg: fitfunc(p_arg, x_arg) - y_arg
    p0 = [0.1,100.0,-3.0]
    p1, success = leastsq(errfunc, p0[:], args=(x, ydata))
    assert success == 1, "Fit did not succeed"
    T2 = 1./p1[1]
    logger.info(strm(T2))
    x_fit = linspace(x.min(),x.max(),5000)
    fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T2 = %0.2f ms)'%(T2*1e3), human_units=False)
    xlabel('t (sec)')
    ylabel('Intensity')
    logger.info(strm("T2:",T2,"s"))
    save_fig = False
    if save_fig:
        savefig('20200108_CPMG_trials.png',
                transparent=True,
                bbox_inches='tight',
                pad_inches=0,
                legend=True)
    fl.show()
