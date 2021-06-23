from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
#from pyspecProcScripts.third_level.process_enhancement import process_enhancement
from sympy import symbols, Symbol, latex,limit,init_printing
#from matplotlib import *
from numpy import *
import matplotlib.pyplot as plt
#from pylab import *
#from sympy import exp as s_exp
#from scipy.optimize import leastsq,minimize,basinhopping,nnls
#import os
#os.chdir('C:/Users/saman/proc_scripts/')
#{{{input parameters
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
fl = fl_mod()
#}}}
s = find_file('210617_T177R1a_pR_DDM_ODNP',exp_type='odnp',expno='enhancement',
        postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=fl)
fl.next('raw data')
s = s['t2':(-1e3,1e3)]
fl.next('sliced to integration bounds')
fl.image(s)
ph0 = s['power',-4].sum('t2')
ph0 /= abs(ph0)
s /= ph0
fl.next('phased')
fl.image(s)
s = s['t2':(-100,100)]
d = s['ph1',1].C
d.integrate('t2')
d /= max(d.data)
#{getting power axis
power_axis_dBm = array(d.get_prop('meter_powers'))
power_axis_W = zeros_like(power_axis_dBm)
power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
power_axis_W = r_[0,power_axis_W]
d.setaxis('power',power_axis_W)
#}
fl.next('simple integral')
fl.plot(d['power',:-3],'ko')
fl.plot(d['power',-3:],'ro')
plt.xlabel('power (W)')
plt.ylabel('integrated 1HNMR signal')
#{Enhancement
E = 1 - d.C
fl.next('enhancement curve')
fl.plot(E['power',:-3],'ko')
fl.plot(E['power',-3:],'ro')
plt.xlabel('power (W)')
fl.show()
#}
