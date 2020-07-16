from pyspecdata import *
from sympy import symbols
def _fitcurve_initial(s,f_range,direct,indirect,guess):
    if not s.get_ft_prop(direct):
        raise ValueError("Your data should be in the frequency domain!")
    rec_curve = s[direct:f_range].sum(direct).real
    print(ndshape(s))
    sgn = sign(rec_curve[indirect:(rec_curve.getaxis(indirect).max())].item())
    rec_curve *= sgn
    f = fitdata(rec_curve)
    return f
def _fitcurve_final(f,whichrate,guess):
    logger.info(strm("Functional form", f.functional_form))
    if guess is not None:
        f.set_guess(guess)
        f.settoguess()
        save_guess = f.eval(100) # if we really wanted to plot the guess,
        # we could return this as well, and then pass it to the fit_curve
        # function
    f.fit()
    print("output:",f.output())
    print("latex:",f.latex())
    if guess is None:
        return f,retval
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
    _fitcurve_initial(s,f_range,direct,indirect,guess)
    M0,Mi,R1,vd = sympy.symbols("M_0 M_inf R_1 %s"%indirect,real=True)
    f.functional_form = Mi + (M0-Mi)*sympy.exp(-vd*R1)
    return _fitcurve_final(f,'R1',guess)

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
    _fitcurve_initial(s,f_range,direct,indirect,guess)
    M0,Mi,R1,vd = sympy.symbols("M_0 R_2 %s"%indirect,real=True)
    f.functional_form = (M0)*sympy.exp(-vd*R1)
    return _fitcurve_final(f,'R1',guess)
