from pyspecdata import *
from pylab import *

with figlist_var() as fl:
    for filename, folder_name, nodename in [
            ('210202_gds_amp_vary_4', 'test_equip', 'capture1')
            ]:
        d = find_file(filename, exp_type=folder_name, expno=nodename)
        fl.next('raw data')
        fl.plot(d)
        d.ft('t',shift=True)
        d = d['t':(0,None)] #toss negative frequencies
        d *= 2 #                multiply data by 2 because the equation
        #                       1/2a*exp(iwt)+aexp(-iwt) and the 2 negated the
        #                       half. taken from analyze_square_refl.py

        fl.next('freq domain')
        fl.plot(d)
        d.ift('t')
        fl.next('analytic signal')
        for j in range(ndshape(d)['amplitudes']):
            fl.plot(abs(d['amplitudes',j]),alpha=0.5, linewidth=1,label = "amplitude%d"%j) 
