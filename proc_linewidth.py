from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
import sympy as s
from symfit import Parameter, Variable, Fit

fl = fl_mod()
for searchstr,exp_type,postproc in [
        ["201118_10mM4AT",'ESR','ESR_linewidth']
        ]:
    d = find_file(searchstr + '.DSC', exp_type=exp_type, postproc=postproc,
                  lookup=postproc_dict)
    fl.next('linewidth for 10mM 4AT')
    fl.plot(d)
    d = d['$B_0$':(-9, 9)]
    fl.next('sliced')
    fl.plot(d)
    s_integral =d.C.run_nopop(cumsum, '$B_0$')
    #fl.next('absorbance')
    #fl.plot(s_integral)
    #{{{fitting with voigt
    B_center = Parameter('B_center')
    sigma = Parameter('sigma')
    gamma = Parameter('gamma')
    A = Parameter('A')
    C = Parameter('C')
    B = Variable('B')
    x = B-B_center
    z = (x+1j*gamma)/(sigma * s.sqrt(2))
    w = (s.exp(-z ** 2)) * s.erfc(-1j * z)
    dVoigt = A * (((gamma * s.im(w)) / ((sigma ** 2) * sigma * s.sqrt(2 * s.pi))) - ((x * s.re(w)) / ((sigma ** 2) * sigma * s.sqrt(2 * s.pi)))) + C
    ydata = d.data.real
    xdata = d.getaxis('$B_0$')
    print(dVoigt)
    #sigma = 0.1
    #gamma = 1.7
    model = dVoigt
    print("MODEL IS")
    print(model)
    fit = Fit(model,xdata,ydata)
    fit_result = fit.execute()
    y = model(B=xdata,A=400,C=fit_result.value(C),sigma=fit_result.value(sigma),gamma=fit_result.value(gamma),B_center=fit_result.value(B_center))
    print(fit_result)
    fl.next('with voigt fit')
    fl.plot(xdata,ydata,'.')
    fl.plot(xdata,y)
    fl.show();quit()
    #{{{fitting to lorentzian
    B0,B,B_center,R,v,v0 = symbols('B0 B B_center R nu nu_0',real=True,positive=True)
    lorentzian = 1/(1j*(B-B_center)+R)
    B_center = Parameter('B_center')
    R = Parameter('R')
    A = Parameter('A')
    #C = Parameter('C')
    B = Variable('B')
    ydata = d.data.real
    xdata = d.getaxis('$B_0$')
    l_deriv = -(A*(2*R*(B+B_center))/(((R**2)+((B+B_center)**2))**2)) 
    print(l_deriv)
    model = l_deriv
    print(model)
    fit = Fit(model,xdata,ydata)
    fit_result = fit.execute()
    y = model(B=xdata,A=fit_result.value(A),B_center=fit_result.value(B_center),
            R=fit_result.value(R))
    print(fit_result.value(R))
    print(fit_result)
    fl.next('fit for 5 mM 4AT')
    fl.plot(xdata,ydata,'.')
    fl.plot(xdata,y)
    #}}}
    fl.show();quit()

