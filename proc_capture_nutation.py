from pyspecdata import *
from pylab import *

with figlist_var() as fl:
    for filename, folder_name, nodename in [
            ('210204_gds_p90_vary_1', 'nutation', 'capture1')
            ]:
        d = find_file(filename, exp_type=folder_name, expno=nodename)
        d.rename('amplitudes','p90')
        print(d['p90',1].get_prop())
        quit()
        fl.next('raw data')
        fl.plot(d)
        d.ft('t',shift=True)
        d = d['t':(0,None)] #toss negative frequencies
        d *= 2 #                multiply data by 2 because the equation
        #                       1/2a*exp(iwt)+aexp(-iwt) and the 2 negated the
        #                       half. taken from analyze_square_refl.py

        fl.next('freq domain')
        fl.plot(d)
        d['t':(None,1.4e7)]=0
        d['t':(1.6e7,None)]=0
        d.ift('t')
        fl.next('analytic signal')
        for j in range(ndshape(d)['p90']):
            fl.plot(abs(d['p90',j]),alpha=0.5, linewidth=1,label = "p90%d"%j)
        fl.next('analytic signal in frequency domain')
        d.ft('t')
        for j in range(ndshape(d)['p90']):
            fl.plot(abs(d['p90',j]),alpha=0.5, linewidth=1,label = "p90%d"%j)

