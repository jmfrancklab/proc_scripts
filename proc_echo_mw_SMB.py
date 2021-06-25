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

s = find_file('210624_Y71R1a_Ras_ODNP',exp_type='ODNP_NMR_comp/ODNP', expno='enhancement',
        postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=fl)
fl.next('raw data')
s = s['t2':(-1e3,1e3)]
fl.image(s)
ph0 = s['power',-4].sum('t2')
ph0 /= abs(ph0)
s /= ph0
fl.next('phased')
fl.image(s)
s = s['t2':(-250,0)]
fl.next('sliced to integration bounds')
fl.image(s)
s = s['ph1',1]
fl.next('sliced to integration bounds\nselect coh. pathway\nreal')
fl.image(s.real)
s.integrate('t2')
s /= max(s.data.real)
#{getting power axis
s.setaxis('power',r_[-9999,array(s.get_prop('meter_powers'))])
print("here are the dBm",s.getaxis('power'))
s.setaxis('power', lambda x:
        1e-3*10**(x/10.))
print("here are the powers",s.getaxis('power'))
s.set_units('power','W')
#}
fl.next('simple integral')
fl.plot(s['power',:-3],'ko')
fl.plot(s['power',-3:],'ro')
plt.ylabel('integrated $^1$H NMR signal')
#{Enhancement
s /= s['power',0]
fl.next('enhancement curve')
fl.plot(s['power',:-3],'ko')
fl.plot(s['power',-3:],'ro')
fl.show()
#}
