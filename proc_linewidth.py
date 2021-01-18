from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
import symfit as s
import pickle,os
from pylab import ndarray
from symfit import Parameter, Variable, Fit
from symfit.core.minimizers import MINPACK
from symfit.contrib.interactive_guess import InteractiveGuess
import numpy as np

fl = fl_mod()
for searchstr,exp_type,postproc in [
        ["210114_50mM_4AT",'ESR','ESR_linewidth']
        ]:
    d = find_file(searchstr + '.DSC', exp_type=exp_type, postproc=postproc,
                  lookup=postproc_dict)
    fl.next('linewidth for 10mM 4AT')
    fl.plot(d)
    d = d['$B_0$':(-10, 10)]
    d.setaxis('$B_0$', lambda x: x+1) # for a positive B_center, b/c the interactive guess doesn't deal well with negative parameters
    fl.next('sliced')
    fl.plot(d)
    s_integral =d.C.run_nopop(np.cumsum, '$B_0$')
    #fl.next('absorbance')
    #fl.plot(s_integral)
    #{{{fitting with voigt
    B_center = Parameter('B_center', value=-0.2)
    sigma = Parameter('sigma', value=3)
    R = Parameter('R', value=7.0)
    A = Parameter('A', value=3e4)
    C = Parameter('C', value=0)
    B = Variable('B')
    y_var = Variable('y')
    if not os.path.exists('dVoigt.pickle'):
        with open('dVoigt.pickle','wb') as fp:
            # cache the expression, which takes some time to generate
            print("no pickle file found -- generating")
            z = ((B-B_center) + s.I*R)/sigma/s.sqrt(2)
            faddeeva = s.simplify(s.exp(-z**2) * s.erfc(-s.I*z))
            voigt = A*s.re(faddeeva)/sigma/s.sqrt(2*s.pi)
            voigt = voigt.simplify()
            dVoigt = voigt.diff(B).simplify()
            pickle.dump(dVoigt,fp)
    else:
        with open('dVoigt.pickle','rb') as fp:
            print("reading expression from pickle")
            dVoigt = pickle.load(fp)
    fl.next('plot guess')
    print(A.value,"a value")
    # {{{ need to re-do b/c defaults are stored in pickle
    B_center.value = 1 
    # setting bounds on B_center fixes the interactive guess, but BFGS chokes on it
    #B_center.min = -5
    #B_center.max = 5
    sigma.value = 0.5
    #sigma.min = 0
    R.value = 1.5
    #R.min = 0
    A.value = 1e3
    #A.min = 0
    C.value = 0
    # }}}
    model_lambda = s.lambdify([B],dVoigt.subs({
        B_center:B_center.value,
        R:R.value,
        A:A.value,
        sigma:sigma.value}),
        modules=[{'ImmutableMatrix': ndarray}, 'numpy', 'scipy'])
    x_finer = r_[d.getaxis('$B_0$')[0]:d.getaxis('$B_0$')[-1]:500j]
    result = model_lambda(x_finer)
    print(type(result),result.shape)
    guess_nddata = nddata(result, [-1], ['$B_0$']).setaxis(
            '$B_0$',x_finer).set_units('$B_0$',d.get_units('$B_0$'))
    fl.plot(d, label='data')
    fl.plot(guess_nddata,':', label='guess')
    fl.next('guess')
    model = s.Model({y_var:dVoigt})
    guess = InteractiveGuess(model, y=d.data.real, B=d.getaxis('$B_0$'), n_points=500)
    guess.execute()
    print(guess)
    y_var = s.Variable('y')
    print("about to run fit")
    fit = s.Fit(model, d.getaxis('$B_0$'), d.data.real)
    fit_result = fit.execute()
    fl.next('data with fit')
    plot(d, '.', label='data')
    fit_nddata = nddata(
            fit.model(B=x_finer, **fit_result.params).y,
            [-1], ['$B_0$']).setaxis('$B_0$', x_finer)
    fl.plot(fit_nddata, label='fit')
    print(fit_result)

fl.show()
