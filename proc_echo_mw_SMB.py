#{imports,etc.
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
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
#}
#{ initializing and symbols
#logger = init_logging("info")
t2 = symbols('t2')
fl = fl_mod()
#}
#{enter parameters and load data
simple_integral = True # don't need to enter t_max or f_range if True.
t_max = 50e-3 # seconds
f_range = (-0.2e3, 0.2e3) # tuple, Hz
#s = find_file('210617_T177R1a_pR_DDM_ODNP',exp_type='odnp',expno='enhancement',
#        postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=fl)
s = find_file('210623_F195R1a_pR_DHPC_ODNP',exp_type='odnp',expno='enhancement',
        postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=None)
#}
#{showing raw data - and cutting off excess
s = s['t2':(-1e3,1e3)]
fl.next('raw data')
fl.image(s)
#}
#{zeroeth order phasing
ph0 = s['power',-4].sum('t2')
ph0 /= abs(ph0)
s /= ph0
fl.next('phased')
fl.image(s)
#}
#{integrating data with bounds
d = s['ph1',1].C
if simple_integral:
    d = d['t2':(-100,100)]
    fl.next('zoomed in')
    fl.image(s['ph1',1][
else:
    d = d['t2',f_range]
    fl.next('sliced to integration bounds')
fl.image(d)
d.integrate('t2')
d /= max(d.data)
#d /= min(d.data)
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
