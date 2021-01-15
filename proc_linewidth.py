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

B_center = Parameter('B_center', value=-0.2)
sigma = Parameter('sigma', value=3)
R = Parameter('R', value=7.0)
A = Parameter('A', value=3e4)
C = Parameter('C', value=0)
B = Variable('B')
y_var = Variable('y')
for searchstr,exp_type,postproc,thisguess,interactive in [
        ("201118_10mM4AT",'francklab_esr/alex','ESR_linewidth',
            { # here, I entered based on the next, and then copied and pasted the result
                A: 952.059548814867,
                B_center: 0.645245462402653,
                R: 1.1452099193249357,
                sigma: 0.7520327837593248,
                },
            False
            ),
        ("201118_1mM4AT",'francklab_esr/alex','ESR_linewidth',
            {

                A:         4e2,
                B_center:  4.220642e-01,
                R:         3.617994e-01,
                sigma:     8.702270e-01,
                },
            False
            )
        ]:
    d = find_file(searchstr + '.DSC', exp_type=exp_type, postproc=postproc,
                  lookup=postproc_dict)
    plt.figure()
    plt.title('linewidth for 10mM 4AT')
    plot(d)
    d = d['$B_0$':(-9, 9)]
    plot(d, '--', alpha=0.5, linewidth=4)
    d.setaxis('$B_0$', lambda x: x+1) # for a positive B_center, b/c the interactive guess doesn't deal well with negative parameters
    s_integral =d.C.run_nopop(np.cumsum, '$B_0$')
    #{{{fitting with voigt
    if not os.path.exists('dVoigt.pickle'):
        with open('dVoigt.pickle','wb') as fp:
            # cache the expression, which takes some time to generate
            print("no pickle file found -- generating")
            z = ((B-B_center) + s.I*R)/sigma/s.sqrt(2)
            faddeeva = s.simplify(s.exp(-z**2) * s.erfc(-s.I*z))
            voigt = A*s.re(faddeeva)/sigma/s.sqrt(2*s.pi)
            voigt *= sigma * R # so adjusting linewidth doesn't change amplitude
            voigt = voigt.simplify()
            # add real below b/c lambdify was giving complex answer
            dVoigt = s.re(s.re(voigt.diff(B)).simplify())
            pickle.dump(dVoigt,fp)
    else:
        with open('dVoigt.pickle','rb') as fp:
            print("reading expression from pickle")
            dVoigt = pickle.load(fp)
    plt.figure()
    plt.title('plot guess')
    print(A.value,"a value")
    # {{{ need to re-do b/c defaults are stored in pickle
    for k,v in thisguess.items():
        k.value = v
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
    plot(d, label='data')
    plot(guess_nddata,':', label='guess')
    model = s.Model({y_var:dVoigt})
    if interactive:
        guess = InteractiveGuess(model, y=d.data.real, B=d.getaxis('$B_0$'), n_points=500)
        guess.execute()
        print(guess)
    y_var = s.Variable('y')
    print("about to run fit")
    fit = s.Fit(model, d.getaxis('$B_0$'), d.data.real)#, minimizer=MINPACK) # really want to use minpack here, but gives "not proper array of floats
    fit_result = fit.execute()
    print("fit is done")
    plt.figure()
    plt.title('data with fit')
    plot(d, '.', label='data')
    fit_nddata = nddata(
            fit.model(B=x_finer, **fit_result.params).y,
            [-1], ['$B_0$']).setaxis('$B_0$', x_finer)
    plot(fit_nddata, label='fit')
    print(fit_result)

plt.show()
