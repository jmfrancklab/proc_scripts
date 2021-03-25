"""This example actually does not work"""
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
B_name = '$B_0$'
thesecolors = cycle(list('bgrcmykw'))
# just not using figure list throughout, to be consistent
#{{{ setting up parameters and parameter lists
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

#}}}
#{{{appending guess lists into parameter list and creating expression list
for j,C in enumerate(C_list):
    A_list.append(Parameter('A%d'%j, value = A_list_guess[j]))
    B_center_list.append(Parameter('B_center%d'%j, value = B_center_list_guess[j]))
    expressions.append(dVoigt.subs({A:A_list[j],B_center:B_center_list[j],
        R:R2+C*k_H,sigma:sigma}))
#}}}
#{{{lambdify the expressions and plot guesses with actual data
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
    plot(datasets[j], label='dataset%d'%j)
    plot(guess_nddata, ':', label='guess%d'%j)
plt.title('guess fit')
#}}}
#{{{defining the model with the expressions
B_list = [Variable('B%d'%j) for j in range(len(C_list))]
y_list = [Variable('y%d'%j) for j in range(len(C_list))]
print("B looks like this:",B)
model = Model({y_list[j]:expressions[j].subs({B:B_list[j]})
    for j in range(len(C_list))})
print([str(B.subs({B:B_list[j]}).atoms(s.Symbol))
    for j in range(len(C_list))])
print([str(expressions[j].subs({B:B_list[j]}).atoms(s.Symbol))
    for j in range(len(C_list))])
kwargs = {str(B_list[j]):datasets[j].getaxis('$B_0$') for j in range(len(C_list))}
kwargs.update(
        {str(y_list[j]):datasets[j].data.real for j in range(len(C_list))}
        )
#}}}
#{{{Fitting the model
fit = Fit(model,
        **kwargs)
use_pickle = True
if not use_pickle:
    fit_result = fit.execute()
    with open('fit_result.pickle','wb') as fp:
        print('generating pickle file')
        pickle.dump(fit_result,fp)
else: 
    assert os.path.exists('fit_result.pickle')
    with open('fit_result.pickle','rb') as fp:
        print('reading fit result from pickle')
        fit_result = pickle.load(fp)
print("fit is done, pickle dumped")
print(fit_result)
x_axes = {'B%d'%j:datasets[j].getaxis('$B_0$') for j in range(len(datasets))}
y_fit = model(
        **fit_result.params,
        **x_axes)
plt.figure()
plt.title('data with fit')
residual_y = []
for j in range(len(datasets)):
    print('y_fit',y_fit[j])
    fit_result = nddata(y_fit[j], [-1], [B_name]
            ).setaxis(B_name, x_axes['B%d'%j])
    thiscolor = next(thesecolors)
    plot(datasets[j],c=thiscolor,alpha=0.5,
            label='data C = %f'%C_list[j])
    plot(fit_result,'--',
            c=thiscolor,alpha=0.5,
            label='fit C = %f'%C_list[j])
    plot(datasets[j]-fit_result,':',c=thiscolor,alpha=0.5,
            label='residual C = %f'%C_list[j])
plt.xlabel('$B_0$/G')
plt.ylabel('Intensity')
plt.legend(**dict(bbox_to_anchor=(1,1),loc=1,borderaxespad=0))
gridandtick(plt.gca())
plt.show()
