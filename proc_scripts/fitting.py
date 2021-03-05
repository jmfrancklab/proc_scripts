from pyspecdata import *
from sympy import symbols
import sympy as sp
import numpy as np
logger=init_logging('info')
def _fitcurve_initial(s,direct):
    if not s.get_ft_prop(direct):
        raise ValueError("Your data should be in the frequency domain!")
    curve = s.real
    return fitdata(curve)
def _fitcurve_final(f,whichrate,guess):
    logger.info(strm("Functional form", f.functional_form))
    logger.info(strm("Functional form", f.functional_form))
    if guess is not None:
        f.set_guess(guess)
        f.settoguess()
        save_guess = f.eval(100) # if we really wanted to plot the guess,
        # we could return this as well, and then pass it to the fit_curve
        # function
    f.fit()
    logger.info(strm("output:",f.output()))
    logger.info(strm("latex:",f.latex()))
    if guess is None:
        return f,1./f.output(whichrate)
    else:
        return f,1./f.output(whichrate),save_guess
def recovery(s,f_range,direct='t2',indirect='indirect',guess=None):
    """Take phased data, slice and integrate it in order to fit an inversion recovery curve

    Parameters
    ----------
    direct: str
        name of the direct dimension
    indirect: str
        name of the indirect dimension, along which the recovery is measured

    Returns
    -------
    f: fitdata
        the full fitdata object with all info about the fit curve
    T1: float
        just the T1 relaxation time
    save_guess: nddata
        returned only if the guess argument is set.
        Give an nddata with all parameters set to the guess value.
    """
    curve = _fitcurve_initial(s,direct)
    sgn = np.sign(curve[indirect:(curve.getaxis(indirect).max())].item())
    curve *= sgn
    M0,Mi,R1,vd = symbols("M_0 M_inf R_1 %s"%indirect,real=True)
    curve.functional_form = Mi + (M0-Mi)*sp.exp(-vd*R1)
    return _fitcurve_final(curve,'R_1',guess)

def decay(s,f_range,direct='t2',indirect='indirect', guess=None):
    """Take phased data, slice and integrate it in order to fit a T2 decay curve

    Parameters
    ----------
    direct: str
        name of the direct dimension
    indirect: str
        name of the indirect dimension, along which the recovery is measured

    Returns
    -------
    f: fitdata
        the full fitdata object with all info about the fit curve
    T1: float
        just the T1 relaxation time
    """
    curve = _fitcurve_initial(s,f_range,direct,indirect,guess)
    sgn = sign(curve[indirect:0].item())
    curve *= sgn
    M0,R2,vd = sympy.symbols("M_0 R_2 %s"%indirect,real=True)
    curve.functional_form = (M0)*sympy.exp(-vd*R2)
    return _fitcurve_final(curve,'R_2',guess)
