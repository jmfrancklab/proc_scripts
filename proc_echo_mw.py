from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
import numpy as np
import matplotlib.pyplot as plt
from sympy import exp as s_exp
rcParams["savefig.transparent"] = True
logger = init_logging("info")
fl = fl_mod()
t2 = symbols('t2')

# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
for searchstr,exp_type,nodename,postproc,freq_range,max_t in [
        #["210311_TEMPOL500uM_DNP_cap_probe_1", 'ODNP_NMR_comp', 'signal',
        #    'spincore_ODNP_v1', (-6000,6000),0.06]
        ["201203_4AT10mM_DNP_cap_probe_1",'ODNP_NMR_comp','signal',
            'spincore_ODNP_v1', (-5000,5000),0.06]
        ]:
    fl.basename = searchstr
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict,fl=fl)
    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range) # visualize the frequency limits
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
    s.ift('t2') # inverse fourier transform into time domain
    logger.debug(strm("THIS IS THE SHAPE"))
    logger.debug(strm(ndshape(s)))
    fl.next('time domain')
    fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    s = slice_FID_from_echo(s,max_t=max_t,fl=None)    # visualize time domain after filtering and phasing
    #{{{apodizing and zero fill
    fl.next('apodize and zero fill')
    R = 5.0/(max_t) # assume full decay by end time
    s *= np.exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1024)
    fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    #}}}
    #{{{select coherence channel in time domain
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s.ft('t2')
    #}}}
    #{{{plotting enhancement curve at lowest and highest powers 
    fl.next('compare highest power to no power')
    idx_maxpower = np.argmax(s.getaxis('power'))
    fl.plot(s['power',0])
    fl.plot(s['power',idx_maxpower])
    #}}}
    #{{{plotting full enhancement curve
    fl.next('full enhancement curve')
    fl.plot(s)
    #}}}
    #{{{plotting enhancement vs power
    fl.next('enhancement')
    enhancement = s['t2':freq_range].sum('t2').real
    enhancement /= enhancement['power',0]
    enhancement.set_units('power','W')
    fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
    fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)
    plt.ylabel('Enhancement')
    #}}}
fl.show();quit()
