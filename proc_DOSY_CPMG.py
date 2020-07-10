from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
import os
init_logging('debug')
fl = figlist_var()
t2 = symbols('t2')
expno = 3
filename = 'ab_jun172020_w0_5'
for searchstr, exp_type, nodename, postproc in [
        ["ab_jun172020_w0_5.*",'test_equip',3,'DOSY_CPMG_v1'],
        ]:
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict)
    s = s['ph8',0]['ph4',1]['m',1]['n',0]
    s = center_CPMG_echo(s)
    s.ft('t2')
    fl.next('abs: request 3')
    fl.image(abs(s))
    fl.next('request 3')
    fl.image(s)
    fl.next('plot indirect 0')
    fl.plot(s['indirect',0])
    fl.show();quit()
    perform_fitting = False
    if perform_fitting:
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
        print("RELATIVE PHASE SHIFT WAS %0.1f\\us and %0.1f$^\circ$", firstorder, angle(zeroorder_rad)/pi*180)
        fl.next('after phased - real ft')
        fl.image(s.real)
        fl.next('after phased - imag ft')
        fl.image(s.imag)
        if s['echo',0].data[:].sum().real < 0:
            s *= -1
        fl.next('plot')
        fl.plot(s)
        print("Ok")
        save_data = False
        if save_data:
            np.savez('proc_DOSY_CPMG_'+filename+'_expno'+str(expno),
                    data = final_spec.data,
                    indirect = final_spec.getaxis('indirect'),
                    echo = final_spec.getaxis('echo'),
                    t2 = final_spec.getaxis('t2'),
                    )
