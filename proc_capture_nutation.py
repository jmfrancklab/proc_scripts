from pyspecdata import *
from pylab import *

with figlist_var() as fl:
    for filename, folder_name, nodename in [
            ('210204_gds_p90_vary_3', 'nutation', 'capture1')
            ]:
        d = find_file(filename, exp_type=folder_name, expno=nodename)
        #print(ndshape(d))
        #quit()
        #d.rename('amplitudes','p90')
        fl.next('raw data')
        fl.plot(d)
        d.ft('t',shift=True)
        d = d['t':(0,None)] #toss negative frequencies
        d *= 2 #                multiply data by 2 because the equation
        #                       1/2a*exp(iwt)+aexp(-iwt) and the 2 negated the
        #                       half. taken from analyze_square_refl.py

        fl.next('freq domain')
        fl.plot(d)
        #fl.show();quit()
        d['t':(None,1.4e7)]=0
        d['t':(1.6e7,None)]=0
        d.ift('t')
        fl.next('analytic signal')
        #took out for loop and hard coding p90 times because only GDS parameters saved over
        # the pp parameters
        fl.plot(abs(d['p90',0]),alpha=0.5, linewidth=1,label = "p90 = 0.1 us")
        fl.plot(abs(d['p90',1]),alpha=0.5,linewidth=1,label = "p90=1.76 us")
        fl.plot(abs(d['p90',2]),alpha=0.5,linewidth=1,label = "p90=3.41 us")
        fl.plot(abs(d['p90',3]),alpha=0.5,linewidth=1,label = "p90=5.07 us")
        fl.plot(abs(d['p90',4]),alpha=0.5,linewidth=1,label = "p90=6.72 us")
        fl.plot(abs(d['p90',5]),alpha=0.5,linewidth=1,label = "p90=8.38 us")
        fl.plot(abs(d['p90',6]),alpha=0.5,linewidth=1,label = "p90=10.03 us")
        fl.plot(abs(d['p90',7]),alpha=0.5,linewidth=1,label = "p90=11.69 us")
        fl.plot(abs(d['p90',8]),alpha=0.5,linewidth=1,label = "p90=13.34 us")
        fl.plot(abs(d['p90',9]),alpha=0.5,linewidth=1,label = "p90=15 us")
        d = abs(d)
        ninety_pulse = d['t':(1.237e-5,3.09e-5)]
        one_eightypulse = d['t':(5.311e-5,8.8e-5)]
        ninety_pulse = ninety_pulse.sum('t')
        fl.next('integrate 90 pulse')
        fl.plot(ninety_pulse,'o')
        one_eightypulse = one_eightypulse.sum('t')
        fl.next('integrate 180 pulse')
        fl.plot(one_eightypulse,'o')

