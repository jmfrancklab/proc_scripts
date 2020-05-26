from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import postproc_dict 
fl = figlist_var()
       
for searchstr, exp_type, nodename, postproc, label_str in [
        #('200221_CPMG_TEMPOLgel_3p0_1','test_equip','signal','CPMG','deadtime=5'),
        #('200221_CPMG_TEMPOLgel_2p9_1','test_equip','signal','CPMG','deadtime=5'),
        #('200304_CPMG_2p6_1','deadtime=5'),
        #('200305_CPMG_3p5_2','deadtime=5'),
        #('200305_CPMG_3p6_2','deadtime=5'),
        ('200305_CPMG_3p7_2','test_equip','signal','CPMG','deadtime=5'),
        #('200305_CPMG_3p7_3','deadtime=5'),
        #('200305_CPMG_3p8_2','deadtime=5'),
        #('200305_CPMG_3p9_2','deadtime=5'),
        #('200305_CPMG_4p0_1','deadtime=5'),
        ]:
    #loads in data and visualizes raw 
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename, postproc=postproc, lookup=postproc_dict)
    nEchoes = s.get_prop('acq_params')['nEchoes']
    fl.next('raw data - chunking ft')
    fl.image(s)
    s.ft(['ph1'])
    fl.next(' image plot coherence-- ft ')
    fl.image(s)
    s.ift('t2')
    
    #moves nScans to be inside other axes
    s.reorder('nScans',first=True)
    
    #select coherence pathway and average nScans
    s = s['ph1',1].C
    s.mean('nScans')

    #move t2 dimension inside
    s.reorder('t2',first=True)

    #centering echo at 0
    echo_center = abs(s)['tE',0].argmax('t2').data.item()
    s.setaxis('t2', lambda x: x-echo_center)
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    fl.next('check center')
    fl.image(s)
    
    #phase correction
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
    
    #minimizing cost function
    sol = basinhopping(costfun, r_[0.,0.],
            minimizer_kwargs={"method":'L-BFGS-B'},
            callback=print_fun,
            stepsize=100.,
            niter=100,
            T=1000.
            )
    zeroorder_rad, firstorder = sol.x
    angle = (zeroorder_rad)/pi*180 
    #applying phase shift
    phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(-1j*2*pi*zeroorder_rad)
    s *= phshift
    logger.info(strm("RELATIVE PHASE SHIFT WAS {:0.1f}\\us and {:0.1f}$^\circ$".format(firstorder,angle)))
    if s['nEchoes',0].data[:].sum().real < 0:
        s *= -1
    logger.info(strm(ndshape(s)))
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)

    #summing along t2 axis and plotting initial decay
    data = s['t2':(-100,150)].sum('t2')
    fl.next('Echo decay')
    x = s.getaxis('nEchoes')
    ydata = data.data.real
    ydata /= max(ydata)
    fl.plot(x,ydata,'-o', alpha=0.7, label='%s'%label_str, human_units=False)
    print(s.dimlabels)
    #fl.show();quit()
    print(ndshape(s))
    print("Beginning T2 curve")
    s_sliced = s['nEchoes',(0,60)]['t2',(-50,100)]
    s = fitdata(s_sliced)
    M0,T2,tE = sympy.symbols("M_0 T_2 t_E", real=True)
    s.functional_form = M0*sympy.exp(-tE/T2)
    print("Functional form", s.functional_form)
    print("Function string",s.function_string)
    s.fit_coeff = r_[-1,1,1]
    fl.next('T2 test')
    fl.plot(s,'o',label=f.name())
    print("symbolic variable:",s.symbolic_vars)
    fl.show();quit()

