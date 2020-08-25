from pyspecdata import *
from sympy import symbols
from proc_scripts import *
import math
from line_profiler import LineProfiler
#@profile
def find_echo_center(s, axis='t2',fl=None):
    """Centers and phases a CPMG echo and returns the centered echo
    
    Parameters
    ----------
    axis: str
        name of the axis you are centering on
        (the direct dimension)

    Returns
    -------
    s: nddata
        contains echo-like data with two or more dimensions
    """
    echo_center = hermitian_function_test(s, fl=fl)
    logger.info(strm("echo center is",echo_center))
    return echo_center 
def center_echo(s, echo_center, axis='t2',fl=None):
    logger.debug(strm("t2 axis start and stop before shifting",s.getaxis(axis)[r_[0,-1]]))
    logger.debug(strm("about to subtract",echo_center,"from the axis"))
    s.setaxis(axis, lambda x: x-echo_center)
    logger.debug(strm("t2 axis start and stop",s.getaxis(axis)[r_[0,-1]]))
    s.register_axis({axis:0})
    s /= zeroth_order_ph(s[axis:0])
    time_bound = min(abs(s.getaxis(axis)[r_[0,-1]]))
    s = s[axis:(-time_bound,time_bound)]
    assert isclose(s.getaxis(axis)[0],-s.getaxis(axis)[-1]),"echo is not symmetric! you are using the wrong code!! (first point is %g, last point %g, dwell time %g, and time_bound=%g"%(s.getaxis(axis)[0],s.getaxis(axis)[-1],diff(s.getaxis(axis)[r_[0,1]]).item(),time_bound)
    return s
def minimize_CPMG_imag(s, axis='t2', fl=None):
    """optimize the first and second order phase of a CPMG pulse sequence
    by minimizing the energy of its imaginary component"""
    #{{{cost function phase correction
    s.ft(axis)
    f_axis = s.fromaxis(axis)
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
