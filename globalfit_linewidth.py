from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
import symfit as s
import pickle,os
from pylab import ndarray
from symfit import Parameter,Variable, Fit, Model
from symfit.core.minimizers import MINPACK
import numpy as np

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

sigma = sympy.symbols('sigma')
sigma = Parameter('sigma', value= 5.9820375e-1)

R2 = sympy.symbols('R2')
R2 = Parameter('R2',value=0.3436)

k_H = sympy.symbols('k_H')
k_H = Parameter('k_H',value=70.302)

R = sympy.symbols('R')
R = Parameter('R', value = 70.6)

A0,A1,A2,A3 = sympy.symbols('A0 A1 A2 A3')
A0 = Parameter('A0' value=2.107721e02)
A1 = Parameter('A1', value = 2.862171e02)
A2 = Parameter('A2', value = 3.527704e02)
A3 = Parameter('A3', value = 4.020849e02)

BC0,BC1,BC2,BC3 = sympy.symbols('BC0 BC1 BC2 BC3')
BC0 = Parameter('BC0', value = 2.391618e-01)
BC1 = Parameter('BC1', value = 1.316801e-01)
BC2 = Parameter('BC2', value = 9.543103e-03)
BC3 = Parameter('BC3', value = 8.16337e-02)

B = sympy.symbols('B')
B = Parameter('B')

for C in enumerate(C_list):
    expression0 = dVoigt.subs({A:A0, B_center:BC0, R: R2 + k_H*C, sigma:sigma})
    expression1 = dVoigt.subs({A:A1, B_center:BC1, R: R2 + k_H*C, sigma:sigma})
    expression2 = dVoigt.subs({A:A2, B_center:BC2, R: R2 + k_H*C, sigma:sigma})
    expression3 = dVoigt.subs({A:A3, B_center:BC3, R: R2 + k_H*C, sigma:sigma})



expressions = []
for C in enumerate(C_list):
    expressions.append(dVoigt.subs({A:A_list_guess[j],B_center:B_center_list_guess[j],
        R:R2+C*k_H,sigma:sigma}))
model_lambda = s.lambdify([B],expressions,
    modules=[{'ImmutableMatrix': ndarray}, 'numpy', 'scipy'])
x_axis = r_[datasets[j].getaxis('$B_0$')[0]:datasets[j].getaxis('$B_0$')[-1]:500j]
print(type(x_axis))
result = model_lambda(x_axis)
#print(type(result),result.shape)
guess_nddata = nddata(result, [-1], ['$B_0$']).setaxis(
        '$B_0$', x_axis).set_units('$B_0$',datasets[j].get_units('$B_0$'))
plot(datasets[j], label='data')
plot(guess_nddata, ':', label='guess')
fl.show();quit()    

        
