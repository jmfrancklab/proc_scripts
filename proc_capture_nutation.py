from pyspecdata import *
from pylab import *

with figlist_var() as fl:
    for filename, folder_name, nodename in [
            ('210202_gds_amp_vary_2', 'nutation', 'capture1')
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
        d['t':(None,1.18e7)]=0
        d['t':(1.8e7,None)]=0
        d.ift('t')
        fl.next('analytic signal')
        #took out for loop and hard coding p90 times because only GDS parameters saved over
        # the pp parameters
        fl.plot(abs(d['amplitudes',0]),alpha=0.5, linewidth=1,label = "p90 = 0 amp")
        fl.plot(abs(d['amplitudes',1]),alpha=0.5,linewidth=1,label = "p90=0.056 amp")
        fl.plot(abs(d['amplitudes',2]),alpha=0.5,linewidth=1,label = "p90=0.111 amp")
        fl.plot(abs(d['amplitudes',3]),alpha=0.5,linewidth=1,label = "p90=0.167 amp")
        fl.plot(abs(d['amplitudes',4]),alpha=0.5,linewidth=1,label = "p90=0.222 amp")
        fl.plot(abs(d['amplitudes',5]),alpha=0.5,linewidth=1,label = "p90=0.278 amp")
        fl.plot(abs(d['amplitudes',6]),alpha=0.5,linewidth=1,label = "p90=0.333 amp")
        fl.plot(abs(d['amplitudes',7]),alpha=0.5,linewidth=1,label = "p90=0.389 amp")
        fl.plot(abs(d['amplitudes',8]),alpha=0.5,linewidth=1,label = "p90=0.444 amp")
        fl.plot(abs(d['amplitudes',9]),alpha=0.5,linewidth=1,label = "p90=0.5 amp")
        fl.next('analytic signal in frequency domain')
        d.ft('t')
        #for j in range(ndshape(d)['p90']):
        #    fl.plot(abs(d['p90',j]),alpha=0.5, linewidth=1,label = "p90%d"%j)

