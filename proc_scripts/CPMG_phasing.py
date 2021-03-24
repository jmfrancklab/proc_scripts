from pyspecdata import *
from sympy import symbols
from proc_scripts import *
import math
import numpy as np
import logging

def center_echo(s, echo_center, axis='t2',fl=None):
    """Slices out a symmetric echo.

    Will generate an error if the echo appears to be extremely lopsided.

    Parameters
    ==========
    echo_center:    float
            center of echo
    axis:           str
            axis along which the echo is being centered
    """        
    s.setaxis(axis, lambda x: x-echo_center)
    s.register_axis({axis:0})
    logging.debug(strm("after register axis, t2 axis is", s.getaxis(axis)))
    s /= zeroth_order_ph(s[axis:0],fl=fl)
    time_bound = min(abs(s.getaxis(axis)[r_[0,-1]]))
    logging.debug(strm("time bound is",time_bound))
    axis_before = ndshape(s)[axis]
    s = s[axis:(-time_bound,time_bound)]
    axis_after = ndshape(s)[axis]
    logging.info(strm(axis_after/axis_before))
    assert axis_after/axis_before > 0.5, "the echo is extremely lopsided -- either you houldn't be using this function, or the center is not actually at echo_center=%g, where you are claiming it is"%echo_center
    logging.info(strm("time bound is",time_bound))
    logging.info(strm("after setting axis to time bounds", s.getaxis(axis)))
    assert np.isclose(s.getaxis(axis)[0],-s.getaxis(axis)[-1]),"echo is not symmetric! you are using the wrong code!! (first point is %g, last point %g, dwell time %g, and time_bound %g"%(s.getaxis(axis)[0],s.getaxis(axis)[-1],np.diff(s.getaxis(axis)[r_[0,1]]).item(),time_bound)
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
        logging.info(strm(iteration, x, f, int(accepted)))
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
    logging.info(strm(ndshape(s)))
    if fl is not None:
        fl.next('after phased - real ft')
        fl.image(s.real)
        fl.next('after phased - imag ft')
        fl.image(s.imag)
    #}}}
    return s
