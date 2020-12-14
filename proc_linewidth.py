from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import *
import sympy as sympy
from symfit import Parameter, Variable, Fit

fl = fl_mod()
for searchstr,exp_type,postproc in [
        ["201214_100mM4AT",'ESR','ESR_linewidth']
        ]:
    s = find_file(searchstr+'.DSC',exp_type=exp_type,postproc=postproc,
            lookup=postproc_dict)
    fl.next('linewidth for 10mM 4AT')
    fl.plot(s)
    s = s['$B_0$':(-11,9)]
    fl.next('sliced')
    fl.plot(s)
    s_integral =s.C.run_nopop(cumsum,'$B_0$')
    #fl.next('absorbance')
    #fl.plot(s_integral)
    B0,B,B_center,R,v,v0 = symbols('B0 B B_center R nu nu_0',real=True,positive=True)
    lorentzian = 1/(1j*(B-B_center)+R)
    B_center = Parameter('B_center')
    R = Parameter('R')
    M = Parameter('M')

    B = Variable('B')
    ydata = s.data.real
    xdata = s.getaxis('$B_0$')
    l_deriv = -M*(2*R*(B+B_center))/(((R**2)+((B+B_center)**2))**2) 
    print(l_deriv)
    model = l_deriv
    print(model)
    fit = Fit(model,xdata,ydata)
    fit_result = fit.execute()
    y = model(B=xdata,M=fit_result.value(M),B_center=fit_result.value(B_center),
            R=fit_result.value(R))
    print(fit_result.value(R))
    print(fit_result)
    fl.next('fit')
    fl.plot(xdata,ydata,'.')
    fl.plot(xdata,y)
    fl.show();quit()

