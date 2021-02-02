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
fl = figlist_var()
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
B_center_list_guess = [2.391618e-01,
        1.316801e-01,
        9.543103e-03,
        8.16337e-02]
B = Parameter('B')
R = Parameter('R')
A = Parameter('A')
B_center = Parameter('B_center')
A_list = []
B_center_list = []
expressions = []
R_expr = []
for j,C in enumerate(C_list):
    A_list.append(Parameter('A%d'%j, value = A_list_guess[j]))
    B_center_list.append(Parameter('B_center%d'%j, value = B_center_list_guess[j]))
    R_expr.append(Parameter('R%d'%j, value = R2 + C*k_H))
    expressions.append(dVoigt.subs({A:A_list[j],B_center:B_center_list[j],
        R:R_expr[j],sigma:sigma}))
print(expressions[0].atoms(sp.Symbol))
for j,C in enumerate(C_list):
    guess_exp_lambda = s.lambdify([B],expressions[j].subs({A:A_list[j].value,
        B_center:B_center_list[j].value,
        R:R_expr[j].value,
        sigma:sigma.value}),
        modules=[{'ImmutableMatrix': ndarray}, 'numpy', 'scipy'])
    x_axis = r_[datasets[j].getaxis('$B_0$')[0]:datasets[j].getaxis('$B_0$')[-1]:500j]
    print(type(x_axis))
    guess = guess_exp_lambda(x_axis)
    print(type(result),result.shape)
    guess_nddata = nddata(result, [-1], ['$B_0$']).setaxis(
            '$B_0$', x_axis).set_units('$B_0$',datasets[j].get_units('$B_0$'))
    plot(datasets[j], label='data')
    plot(guess_nddata, ':', label='guess')
fl.show();quit()    

        
