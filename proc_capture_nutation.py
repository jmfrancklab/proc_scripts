from pyspecdata import *
from pylab import *
import matplotlib.pyplot as plt
with figlist_var() as fl:
    for filename, folder_name, nodename in [
            ('210204_gds_p90_vary_3', 'nutation', 'capture1')
            ]:
        d = find_file(filename, exp_type=folder_name, expno=nodename)
        #print(ndshape(d))
        #quit()
        #d.rename('amplitudes','p90')
        #d = d['t':(1.2e-5,3e-5)]
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
        d['t':(None,1.26e7)]=0
        d['t':(1.73e7,None)]=0
        d.ift('t')
        #fl.next('analytic signal')
        #took out for loop and hard coding p90 times because only GDS parameters saved over
        # the pp parameters
        d = d.sum('t')
        d0 = np.array(abs(d['p90':0]))
        d1 = abs(d['p90',1])
        d2 = abs(d['p90',2])
        d3 = abs(d['p90',3])
        d4 = abs(d['p90',4])
        d5 = abs(d['p90',5])
        d6 = abs(d['p90',6])
        d7 = abs(d['p90',7])
        d8 = abs(d['p90',8])
        d9 = abs(d['p90',9])
        print(d3)
        print(d4)
        print(d5)
        print(d6)
        print(d7)
        print(d8)
        print(d9)
        plt.figure()
        plt.title('tip angles')
        plt.plot(0,4.901562377307104e-19,'ro')
        plt.plot(1.67e-6,3.6088313067019186e-17,'ro')
        plt.plot(3.33e-6,4.758550126787486e-17,'ro')
        plt.plot(5e-6, 3.957033682206997e-17,'ro')
        plt.plot(6.67e-6,1.1355484503437166e-16,'ro')
        plt.plot(8.33e-6,4.5102810375396984e-17,'ro')
        plt.plot(10e-6,7.006557745786588e-17,'ro')
        plt.plot(11.67e-6,1.2213488147284855e-16,'ro')
        plt.plot(13.33e-6,1.2860667056498242e-16,'ro')
        plt.plot(15e-6,1.090429294243217e-16,'ro')
        plt.xlabel('time')
        plt.ylabel('abs(d[p90,j]).sum(t)')
        show()

        #fl.plot(abs(d)['p90':0],alpha=0.5, linewidth=1,label = "p90 = 0 us")
        #fl.plot(abs(d['p90',1]),alpha=0.5,linewidth=1,label = "p90= 1.67 us")
        #fl.plot(abs(d['p90',2]),alpha=0.5,linewidth=1,label = "p90= 3.33 us")
        #fl.plot(abs(d['p90',3]),alpha=0.5,linewidth=1,label = "p90= 5 us")
        #fl.plot(abs(d['p90',4]),alpha=0.5,linewidth=1,label = "p90= 6.67 us")
        #fl.plot(abs(d['p90',5]),alpha=0.5,linewidth=1,label = "p90= 8.33 us")
        #fl.plot(abs(d['p90',6]),alpha=0.5,linewidth=1,label = "p90= 10 us")
        #fl.plot(abs(d['p90',7]),alpha=0.5,linewidth=1,label = "p90= 11.67 us")
        #fl.plot(abs(d['p90',8]),alpha=0.5,linewidth=1,label = "p90= 13.33 us")
        #fl.plot(abs(d['p90',9]),alpha=0.5,linewidth=1,label = "p90= 15 us")
        #fl.next('analytic signal in frequency domain')
        #for j in range(len(ndshape(d)['p90'])):
        #    fl.plot(abs(d['p90',j]),alpha=0.5, linewidth=1,label = "p90%d"%j)

