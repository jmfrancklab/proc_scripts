from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.correlation_alignment_ODNP import correl_align
from sympy import symbols
from matplotlib import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from sympy import exp as s_exp
from itertools import cycle
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
thesecolors = cycle(list('bgrcmykw'))
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
def as_scan_nbr(s):
        return s.C.setaxis('power','#').set_units('power','scan #')
# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
def process_enhancement(s, searchstr='', signal_pathway = {'ph1':1,'ph2':-2},
        excluded_pathways = [(0,1),(0,0)], freq_range=(None,None),
        t_range=(0,0.083),fl=None):
    if fl is not None:
        fl.basename = searchstr
        fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
                s,freq_range) # visualize the frequency limits
    s.ift('t2')
    if fl is not None:
        fl.basename = searchstr
        fl.side_by_side('show time limits',s,t_range)
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    if fl is not None:
        fl.next('look at data')
        fl.image(s.C.setaxis('power','#'))
    s.ift(['ph1','ph2'])    
    rx_offset_corr = s['t2':(0.008,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('power',raw_index=True).item()
    s = s['t2':freq_range]    
    s.ft(['ph1','ph2'])
    if fl is not None:
        fl.next('')
        fl.image(s.C.setaxis('power','#'))
    s.ift('t2') # inverse fourier transform into time domain
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    logger.info(strm("applying zeroth order correction"))
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',1]['ph2',0] = 1
    phasing.ift(['ph1','ph2'])
    if fl is not None:
        fl.next('zeroth order corrected')
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
    s /= ph0
    s.ft(['ph1','ph2'])
   #{{{apodizing and zero fill
    logger.info(strm(s.dimlabels))
    if fl is not None:
        fl.next('phase corrected')
        fl.image(as_scan_nbr(s))
    s.reorder(['ph1','ph2','power','t2'])
    logger.info(strm("zero corssing at",zero_crossing))
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = 1e-3*10**(power_axis_dBm/10)
    power_axis_W = r_[0,power_axis_W]
    s.ift(['ph1','ph2'])
    s.ft('t2')
    frq_max = abs(s).argmax('t2')
    s.ift('t2')
    s *= np.exp(-1j*2*pi*frq_max*s.fromaxis('t2'))
    s.ft('t2')
    if fl is not None:
        fl.next(r'after rough alignment, $\varphi$ domain')
        fl.image(as_scan_nbr(s))
    fl.basename='correlation subroutine --before zero crossing:'
    logger.info(strm("ndshape",ndshape(s),"zero crossing at",zero_crossing))
    if zero_crossing > 1:
        opt_shift,sigma = correl_align(s['power',:zero_crossing+1],indirect_dim='power',
                ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
                sigma=50)
        s.ift('t2')
        s['power',:zero_crossing+1] *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
        s.ft('t2')
    else:
        logger.warning("You have 1 point or less before your zero crossing!!!!")
    fl.basename = 'correlation subroutine -- after zero crossing'
    opt_shift,sigma = correl_align(s,indirect_dim='power',
            ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
            sigma=100)
    s.ift('t2')
    s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename = None
    if fl is not None:
        fl.next(r'after correlation, $\varphi$ domain')
        fl.image(as_scan_nbr(s))
    s.ift('t2')
    s.ft(['ph1','ph2'])
    if fl is not None:
        fl.next('after correlation alignment FTed ph')
        fl.image(as_scan_nbr(s))
    s.reorder(['ph2','ph1','power','t2'])
    if fl is not None:
        fl.next('after correlation -- time domain')
        fl.image(as_scan_nbr(s))
    s.ft('t2')    
    if fl is not None:
        fl.next('after correlation -- frequency domain')
        fl.image(as_scan_nbr(s))
    s.ift('t2')
    s = s['t2':(0,None)]
    s['t2':0] *= 0.5
    s.ft('t2')
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    error_pathway = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    # }}}
    s_,frq_slice = integral_w_errors(s,signal_pathway,error_pathway,
            indirect='power', fl=fl, return_frq_slice=True)
    s = s_.C
    idx_maxpower = np.argmax(s.getaxis('power'))
    if fl is not None:
        fl.next('full enhancement curve')
        fl.plot(s)
    s /= max(s.data)
    s.setaxis('power',power_axis_W)
    thiscolor = next(thesecolors)
    if fl is not None:
        fl.next('E(p)')
        fl.plot((s['power',:idx_maxpower+1]), 'ko', capsize=6,
                alpha=0.3)
        fl.plot((s['power',idx_maxpower+1:]),'ro',capsize=6,
                alpha=0.3)
    enhancement = s
    show()
    return enhancement,idx_maxpower
