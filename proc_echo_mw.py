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
freq_range= (-5e3,5e3)
t_range = (0,50e-3)
fl = fl_mod()
#}}}
for filename,nodename,file_location,postproc in [
        ('210615_S175R1a_pR_DDM_ODNP','enhancement',
            'odnp','spincore_ODNP_v1'),
        ]:
    s = find_file(filename,exp_type=file_location,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    myslice = s['t2':freq_range]
    mysign = determine_sign(select_pathway(myslice,signal_pathway,mult_ph_dims=False))
    enhancement,idx_maxpower = process_enhancement(s,searchstr = filename,
            freq_range=freq_range, t_range=t_range,flip=True,sign=mysign,fl=fl)
    fl.show()

