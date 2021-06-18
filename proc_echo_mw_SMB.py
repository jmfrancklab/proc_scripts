from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from pyspecProcScripts.third_level.process_enhancement import process_enhancement
from sympy import symbols, Symbol, latex,limit,init_printing
#from matplotlib import *
#import numpy as np
#import matplotlib.pyplot as plt
#from pylab import *
#from sympy import exp as s_exp
#from scipy.optimize import leastsq,minimize,basinhopping,nnls
import os
os.chdir('C:/Users/saman/proc_scripts/')
#{{{input parameters
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
signal_pathway = {'ph1':1}
fl = fl_mod()
#}}}
s = find_file('210617_T177R1a_pR_DDM_ODNP',exp_type='odnp',expno='enhancement',
        postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=fl)
print(s.get_prop('acq_params'));quit()
fl.next('raw data')
s = s['t2':(-1e3,1e3)]
fl.next('sliced to integration bounds')
fl.image(s)
ph0 = s['power',-4].sum('t2')
ph0 /= abs(ph0)
s /= ph0
fl.image(s)
fl.next('simple integral')
s = s['t2':(-100,100)]
s.integrate('t2')
fl.plot(s,'o')
fl.show()
print(s.getaxis('power'))
