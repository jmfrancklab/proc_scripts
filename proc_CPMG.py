from pyspecdata import *
from scipy.optimize import basinhopping
from proc_scripts import postproc_dict 
fl = figlist_var()
logger = init_logging('info')

for searchstr, exp_type, nodename, postproc, label_str in [
        #('200221_CPMG_TEMPOLgel_3p0_1','test_equip','signal','CPMG','deadtime=5'),
        ('200221_CPMG_TEMPOLgel_2p9_1','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200304_CPMG_2p6_1','test_equip','signal','CPMG','deadtime=5'),
        #('200305_CPMG_3p5_2','deadtime=5'),
        #('200305_CPMG_3p6_2','deadtime=5'),
        #('200305_CPMG_3p7_2','test_equip','signal','CPMG','deadtime=5'),
        #('200305_CPMG_3p7_3','deadtime=5'),
        #('200305_CPMG_3p8_2','deadtime=5'),
        #('200305_CPMG_3p9_2','deadtime=5'),
        #('200305_CPMG_4p0_1','deadtime=5'),
        ]:
    ###{{{loading in data and displaying raw data
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename, postproc=postproc, lookup=postproc_dict)
    #}}}
    #{{{select and display coherence channel centered
    s.ft(['ph1'])
    s.ift('t2')
    s.reorder('nScans',first=True)
    s = s['ph1',1].C
    s.mean('nScans')
    s.reorder('t2',first=True)
    echo_center = abs(s)['tE',0].argmax('t2').data.item()
    s.setaxis('t2', lambda x: x-echo_center)
    fl.next('check center')
    fl.image(s)
    #}}}
    #{{{cost function phase correction
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
        logger.info(strm(iteration, x, f, int(accepted)))
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
    logging.info(strm("RELATIVE PHASE SHIFT WAS %0.1f\\us and %0.1f$^\circ$", firstorder, angle(zeroorder_rad)/pi*180))
    if s['tE',0].data[:].sum().real < 0:
        s *= -1
    logger.info(strm(ndshape(s)))
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)
    #}}}
    #{{{select echo decay fit function
    s.ift('t2')
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    data = s['t2':(0,None)].sum('t2')
    fl.next('Echo decay')
    fl.plot(data,'o')
    print("starting T2 curve")
    f = fitdata(data)
    M0,R2,tE = sympy.symbols("M_0 R_2 tE", real=True)
    f.functional_form = M0*sympy.exp(-tE*R2)
    fl.next('T2 test')
    fl.plot(f,'o',label=f.name())
    f.fit()
    fl.plot(f.eval(100),label='%s fit'%f.name())
    text(0.75,0.25, f.latex(),transform=gca().transAxes, size='large',
            horizontalalignment='center', color= 'k')
    print("output",f.output())
    print("latex",f.latex())
    T2 = 1./f.output('R_2')
    fl.show();quit()
    #}}}
    #{{{saving figure
    save_fig = False
    if save_fig:
        savefig('20200108_CPMG_trials.png',
                transparent=True,
                bbox_inches='tight',
                pad_inches=0,
                legend=True)
    fl.show()
    #}}}
