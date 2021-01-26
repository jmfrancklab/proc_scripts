from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
import symfit as s
import pickle,os
from pylab import ndarray
from symfit import Parameter, Variable, parameters, variables, Fit, Model
from symfit.core.minimizers import MINPACK
from symfit.contrib.interactive_guess import InteractiveGuess
import numpy as np
from itertools import cycle
from scipy.optimize import nnls

B_center = Parameter('B_center', value=-0.2)
sigma = Parameter('sigma', value=3)
R = Parameter('R', value=7.0)
A = Parameter('A', value=3e4)
C = Parameter('C', value=0)
B = Variable('B')
y_var = Variable('y')
conc_list = []
R_list = []
sigma_list = []
data = []
A_list = []
B_center_list = []
C_list = []
for searchstr,exp_type,postproc,thisguess,interactive,concentration in [
        ("210114_3mM_4AT",'ESR','ESR_linewidth',
            {
                A:2.107721e02,
                B_center: 2.391618e-01,
                R: 5.076452e-01,
                sigma: 6.672132e-01,
                },
            False,
            3e-3
            ),
        ("210114_5mM_4AT",'ESR','ESR_linewidth',
            {
                A:2.862171e02,
                B_center: 1.316801e-01,
                R: 7.333833e-01,
                sigma: 5.831063e-01,
                },
            False,
            5e-3
            ),
        ("210114_7mM_4AT",'ESR','ESR_linewidth',
            {
                A: 3.527704e02,
                B_center: 9.543103e-03,
                R: 8.813178e-01,
                sigma: 5.54451e-01,
                },
            False,
            7e-3
            ),
        ("210114_10mM_4AT",'ESR','ESR_linewidth',
            { # here, I entered based on the next, and then copied and pasted the result
                A: 4.020849e02,
                B_center: 8.16337e-02,
                R: 1.009644,
                sigma: 5.880406e-01,
                },
            False,
            10e-3
            ),
        ("201118_1mM4AT",'ESR','ESR_linewidth',
            {

                A:         4e2,
                B_center:  4.220642e-01,
                R:         3.617994e-01,
                sigma:     8.702270e-01,
                },
            False,
            1e-3
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
    conc_list.append(concentration)
    R_list.append(fit_result.params['R'])
    sigma_list.append(fit_result.params['sigma'])
    A_list.append(fit_result.params['A'])
    B_center_list.append(fit_result.params['B_center'])
    data.append(d)
plt.figure()
print(data)
print(data[0].dimlabels)
plot(conc_list,R_list,'x',label='R')
plot(conc_list,sigma_list,'x',label=r'$\sigma$')
plt.legend(**dict(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.))
plt.ylabel('concentration')
plt.show()

#{{{ attempting global fitting
thesecolors = cycle(list('bgrcmykw'))
A = parameters('A')
R = parameters('R')
sigma = parameters('sigma')
#B = parameters('B')
B_center = parameters('B_center')
y0,y1,y2,y3 = variables('y0, y1, y2, y3')
B0_0,B0_1,B0_2,B0_3 = variables('B0_0, B0_1, B0_2, B0_3')
#C = variables('C0, C1, C2, C3')
model = Model({
    y0: dVoigt.subs({A:A_list[0]}).subs({R:R_list[0]}).subs({sigma:sigma_list[0]}).subs({B_center:B_center_list[0]}),
    y1: dVoigt.subs({A:A_list[1]}).subs({R:R_list[1]}).subs({sigma:sigma_list[1]}).subs({B_center:B_center_list[1]}),
    y2: dVoigt.subs({A:A_list[2]}).subs({R:R_list[2]}).subs({sigma:sigma_list[2]}).subs({B_center:B_center_list[2]}),
    y3: dVoigt.subs({A:A_list[3]}).subs({R:R_list[3]}).subs({sigma:sigma_list[3]}).subs({B_center:B_center_list[3]}),
    })
fit = Fit(model,
        y0 = data[0].data.real,
        y1 = data[1].data.real,
        y2 = data[2].data.real,
        y3 = data[3].data.real,
        #B0_0 = data[0].getaxis('$B_0$'),
        #B0_1 = data[1].getaxis('$B_0$'),
        #B0_2 = data[2].getaxis('$B_0$'),
        #B0_3 = data[3].getaxis('$B_0$'),
        )
fit_result = fit.execute()
x_finer = r_[0:10.0:500j]
y_fit = model(
        A = fit_result.value(A),
        R = fit_result.value(R),
        sigma = fit_result.value(sigma),
        B_center = fit_result.value(B_center),
        B0_0 = x_finer,
        B0_1 = x_finer,
        B0_2 = x_finer,
        B0_3 = x_finer,
        )
A = fit_result.value(A)
R = fit_result.value(R)
sigma = fit_result.value(sigma)
B_center = fit_result.value(B_center)

fl.next('global fit')
for j in range(3):
    thiscolor = next(thesecolors)
    fl.plot(data[j].getaxis('$B_0$'),data[j].real,'.',c=thiscolor,alpha=0.7)
    fl.plot(x_finer,y_fit[j],':',c=thiscolor,alpha=0.7,label='symfit,global')
gridandtick(gca())
legend()
fl.show();quit()
