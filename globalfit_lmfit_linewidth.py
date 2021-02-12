from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
import symfit as s
import pickle,os
from pylab import ndarray
from symfit import Parameter,Variable, Fit, Model
from symfit.core.minimizers import MINPACK
import numpy as np
import sympy as sp
from itertools import cycle
from lmfit import Parameters, minimize,Minimizer
B_name = '$B_0$'
thesecolors = cycle(list('bgrcmykw'))
fl = figlist_var()
#{{{ starting with sympy expression
datasets = []
C_list = []
for thisfile,C in [
        ('210114_3mM_4AT', 3e-3),
        ('210114_5mM_4AT', 5e-3),
        ('210114_7mM_4AT', 7e-3),
        ('210114_10mM_4AT', 10e-3)
        ]:
    datasets.append(find_file(thisfile +'.DSC', exp_type='ESR', postproc='ESR_linewidth',
        lookup=postproc_dict))
    C_list.append(C)
with open('dVoigt.pickle','rb') as fp:
    print("reading expression for dVoigt from pickle")
    dVoigt = pickle.load(fp)
sigma = Parameter('sigma', value= 5.9820375e-1)
R2 = Parameter('R2', value = 0.3436)
k_H = Parameter('k_H', value = 70.302)
A_list_guess = [2.107721e02,
        2.862171e02,
        3.527704e02,
        4.020849e02]
B_center_list_guess = [2.391618e-01-1,# the -1 here comes from the fact that we mess with the x axis in the proc_linewidth code where these come from
        1.316801e-01-1,
        9.543103e-03-1,
        8.16337e-02-1]
B = Variable('B')
R = Parameter('R')
A = Parameter('A')
B_center = Parameter('B_center')
A_list = []
B_center_list = []
expressions = []

for j,C in enumerate(C_list):
    A_list.append(Parameter('A%d'%j, value = A_list_guess[j]))
    B_center_list.append(Parameter('B_center%d'%j, value = B_center_list_guess[j]))
    expressions.append(dVoigt.subs({A:A_list[j],B_center:B_center_list[j],
        R:R2+C*k_H,sigma:sigma}))
for j,C in enumerate(C_list):
    guess_exp_lambda = s.lambdify([B],expressions[j].subs({R2:R2.value,
        A_list[j]:A_list[j].value,
        k_H:k_H.value,
        sigma:sigma.value,
        B_center_list[j]:B_center_list[j].value}),
        modules=[{'ImmutableMatrix': ndarray}, 'numpy', 'scipy'])
    x_axis = r_[datasets[j].getaxis('$B_0$')[0]:datasets[j].getaxis('$B_0$')[-1]:5280j]
    print(type(x_axis))
    guess = guess_exp_lambda(x_axis)
    print(type(guess),guess.shape)
    guess_nddata = nddata(guess, [-1], ['$B_0$']).setaxis(
            '$B_0$', x_axis).set_units('$B_0$',datasets[j].get_units('$B_0$'))
    #}}}
#{{{starting lmfit attempt
p_true = Parameters()
B = datasets[j].getaxis('$B_0$')

for j,C in enumerate(C_list):
    p_true.add('A%d'%j, value = A_list_guess[j])
    p_true.add('B_center%d'%j,value = B_center_list_guess[j])
    p_true.add('sigma',value=5.9820375e-1)
    p_true.add('R2', value= 0.3436)
    p_true.add('k_H',value=70.302)
    p_true.add('R')
    def residual(pars, B, data=None):
        model = guess_nddata
        if data is None:
            return model
        return model - datasets[j]
fit_params = Parameters()
print(expressions[1].atoms(s.Symbol))
for j,C in enumerate(C_list):
    fit_params.add('A%d'%j)
    fit_params.add('B_center%d'%j)
    fit_params.add('sigma')
    fit_params.add('R2')
    fit_params.add('k_H')
    fit_params.add('B')
    out = minimize(residual, fit_params, args=(B,), kws={'data':datasets[j]})
