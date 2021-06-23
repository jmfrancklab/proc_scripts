"""
Processes enhancement data
==========================
Processes data acquired from an enhancement experiment 
and plots the resulting enhancement curve normalized.
"""
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from sympy import symbols
from matplotlib import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from sympy import exp as s_exp
from itertools import cycle
from pyspecProcScripts.third_level.process_data import proc_data
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
thesecolors = cycle(list('bgrcmykw'))
def as_scan_nbr(s):
        return s.C.setaxis(
'nScans','#').set_units('nScans','scan #').setaxis(
'power','#').set_units('power','scan #')
# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
def process_enhancement(s, signal_pathway = {'ph1':1},
        excluded_pathways = [(0,0)], freq_range=(None,None),
        t_range=(0,0.083),direct='t2',sign=None,fl=None):
    """
    Parameters
    ==========
    s:  nddata
    signal_pathway: dict
                    Coherence transfer pathway in which the signal
                    resides.
    excluded_pathways:  lst
                        The coherence transfer pathway in which we
                        do NOT expect to see signal
    freq_range: tuple
                Range in the frequency domain in which the signal resides.
                Units should be Hz.
    t_range:    tuple
                Range in the time domain in which the signal lasts.
                Units should be seconds.
    direct:     str
                Direct dimension.
    """            
                    
    if fl is not None:
        fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
                s,thisrange=freq_range) # visualize the frequency limits
    s.ift('t2')
    s.reorder(['ph1','power','t2'])
    if fl is not None:
        fl.push_marker()
        fl.next('time domain')
        fl.image(s)
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s = proc_data(s,indirect='power',fl=fl, 
            signal_pathway = signal_pathway,f_range=freq_range,
            t_range=t_range,sign=sign)
    #{{{Normalizing by max
    d = s.C
    d /= max(d.data)
    #}}}
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    #power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
    #power_axis_W = r_[0,power_axis_W]
    #d.setaxis('power',power_axis_W)
    thiscolor = next(thesecolors)
    if fl is not None:
        fl.next('E(p)')
        fl.plot(d['power',:-3], 'ko', capsize=6, alpha=0.3)
        fl.plot(d['power',-3:],'ro',capsize=6, alpha=0.3)
        fl.pop_marker()
    enhancement = d
    return enhancement
