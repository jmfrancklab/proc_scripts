from pylab import *
from pyspecdata import *
from scipy.optimize import minimize, leastsq
from sympy import exp as s_exp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, latex, Symbol
from pyspecProcScripts import *
from .simple_functions import select_pathway
from .correlation_alignment import correl_align
from .integral_w_error import integral_w_errors
from .DCCT_func import DCCT
this_figsize = (6,12)
t2 = symbols('t2')
logger = init_logging('info')
def process_data(s,searchstr='',
        signal_pathway={'ph1':0,'ph2':1},
        excluded_pathways = [(0,0)],
        f_range=(None,None),
        t_range=(0,0.083),
        sgn=None,
        direct='t2',
        indirect='indirect',
        error_bars = True,
        correlate = True,
        fl=None):
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    s *= sgn
    s.ift(direct)
    if indirect is 'vd':
        s['ph2',0]['ph1',0]['t2':0] = 0 #kill axial noise
        s.reorder(['ph1','ph2','vd','t2'])
    else:
        s['ph1',0]['t2':0] = 0
        s.reorder(['ph1','power','t2'])
    #{{{DC offset correction
    s.ift(list(signal_pathway))
    t_start = t_range[-1] / 4
    t_start *= 3
    rx_offset_corr = s['t2':(t_start,None)]
    rx_offset_corr = rx_offset_corr.data.mean()
    s -= rx_offset_corr
    s.ft(direct)
    s.ft(list(signal_pathway))
    #}}}
    #{{{phase correction
    s = s[direct:f_range]
    s.ift(list(signal_pathway))
    if fl is not None:
        fl.basename = "(%s)"%searchstr
        DCCT(s,fl.next('%s'%searchstr,figsize=this_figsize),
                total_spacing=0.2)
    s.ft(list(signal_pathway))
    s.ift(direct)
    best_shift = hermitian_function_test(select_pathway(s.C.mean(indirect),signal_pathway).C.convolve(direct,3e-4))
    logger.info(strm("best shift is", best_shift))
    s.setaxis(direct,lambda x: x-best_shift).register_axis({direct:0})
    s /= zeroth_order_ph(select_pathway(s,signal_pathway))
    if fl is not None:
        s.ft(direct)
        DCCT(s,fl.next('Phased Data for %s'%searchstr,figsize=this_figsize),
                total_spacing=0.2)
        s.ift(direct)
    if indirect is 'vd':
        s.reorder(['ph1','ph2','vd',direct])
    else:
        s.reorder(['ph1',indirect,direct])
    #}}}
    #{{{Correlate
    if correlate:
        s.ft(direct)
        s.ift(list(signal_pathway))
        fl.basename='correlation subroutine:'
        opt_shift, sigma, my_mask = correl_align(s,indirect_dim=indirect,
                signal_pathway=signal_pathway,sigma=125)
        s.ift(direct)
        s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis(direct))
        s.ft(direct)
        fl.basename=None
        s.ift(direct)
        s.ft(list(signal_pathway))
        s.ft(direct)
        if indirect is 'vd':
            s.reorder(['ph1','ph2',indirect,direct])
        else:
            s.reorder(['ph1',indirect,direct])
        if fl is not None:
            DCCT(s,fl.next('Aligned %s'%searchstr,figsize=this_figsize),
                    total_spacing=0.2)
     #}}}    
    s.ift(direct)
    s = s[direct:(0,t_range[-1])]
    s[direct,0] *= 0.5
    s.ft(direct)
    s *= sgn
    if fl is not None:
        DCCT(s,fl.next('FID sliced %s'%searchstr,figsize=this_figsize),
                total_spacing=0.2)
    #{{{Integrate with error bars
    if error_bars:
        if indirect is 'vd':
            error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
                    - set(excluded_pathways)
                    - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
            error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
        else:
            error_path = (set(((j) for j in range(ndshape(s)['ph1'])))
                    - set(excluded_pathways)
                    - set([(signal_pathway['ph1'])]))
            error_path = [{'ph1':j} for j in error_path]
        s_int,frq_slice = integral_w_errors(s,signal_pathway,error_path,
                indirect=indirect, return_frq_slice = True)
        x = s_int.get_error()
        x[:] /= sqrt(2)
        if indirect is 'vd':
            s_int = s_int
        else:
            idx_maxpower = np.argmax(s.getaxis('power'))
            s_int *= -1
            s_int /= max(s_int.data)
    #}}}    
    else:
        s_int = s.mean(direct)
        s_int = select_pathway(s_int,signal_pathway)
    if fl is not None:
        fl.basename=None
        fl.next('Integrated Data for %s'%searchstr)
        fl.plot(s_int,'o')
    return s_int, s
    
