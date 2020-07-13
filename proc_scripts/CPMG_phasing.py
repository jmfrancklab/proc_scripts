from pyspecdata import *
from sympy import symbols
from proc_scripts import *
import math
fl=figlist_var()
def center_CPMG_echo(s):
    echo_center = hermitian_function_test(s,fl=fl)#['indirect',0],fl=fl)
    logger.info("echo center is",echo_center)
    s.setaxis('t2', lambda x: x-echo_center)
    s.register_axis({'t2':0})
    s /= zeroth_order_ph(s['t2':0])
    time_bound = min(abs(s.getaxis('t2')[r_[0,-1]]))
    s = s['t2':(-time_bound,time_bound)]
    assert math.isclose(s.getaxis('t2')[0],-s.getaxis('t2')[-1], rel_tol=0.05),"echo is not symmetric! you are using the wrong code!!"
    return s
def minimize_CPMG_imag(s, fl=None):
    """optimize the first and second order phase of a CPMG pulse sequence
    by minimizing the energy of its imaginary component"""
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
    if fl is not None:
        fl.next('after phased - real ft')
        fl.image(s.real)
        fl.next('after phased - imag ft')
        fl.image(s.imag)
    #}}}
    return s
