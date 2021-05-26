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
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
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
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
    s.ift('t2')
    if fl is not None:
        fl.next('look at data')
        fl.image(s.C.setaxis('power','#'))
    rx_offset_corr = s['t2':(0.02,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    if fl is not None:
        fl.next('')
        fl.image(s.C.setaxis('power','#'))
    s.ift('t2') # inverse fourier transform into time domain
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    coh_slice = select_pathway(s['t2':0],signal_pathway)
    s.ft('t2')
    if fl is not None:
        fl.next('before zeroth order phasing correction')
        fl.image(as_scan_nbr(s))
    s.ift('t2')
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
    s.ft('t2')
    if fl is not None:
        fl.next('phase corrected')
        fl.image(as_scan_nbr(s))
    s.reorder(['ph1','ph2','power','t2'])
    zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('power',raw_index=True).item()
    logger.info(strm("zero corssing at",zero_crossing))
    s.ift(['ph1','ph2'])
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
    opt_shift,sigma = correl_align(s['power',zero_crossing+1:],indirect_dim='power',
            ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
            sigma=50)
    s.ift('t2')
    s['power',zero_crossing+1:] *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename = None
    if fl is not None:
        fl.next(r'after correlation, $\varphi$ domain')
        fl.image(as_scan_nbr(s))
    s.ft(['ph1','ph2'])
    if fl is not None:
        fl.next('after correlation alignment FTed ph')
        fl.image(as_scan_nbr(s))
    s.reorder(['ph2','ph1','power','t2'])
    if fl is not None:
        fl.next('after correlation -- frequency domain')
        fl.image(as_scan_nbr(s))
    s.ift('t2')
    if fl is not None:
        fl.next('after correlation -- time domain')
        fl.image(as_scan_nbr(s))
    s = s['t2':(0,t_range[-1])]
    s['t2':0] *= 0.5
    if select_pathway(s['t2':0]['power',0],signal_pathway).real < 0:
        s *= -1
    R = 5.0/(t_range[-1]) # assume full decay by end time
    s *= np.exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=2048)
    if fl is not None:
        fl.next('apodized and zero fill')
        fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    #}}}
    #{{{select coherence channel in time domain
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
    # }}}
    enhancement = integral_w_errors(s,signal_pathway,error_path,indirect='power')
     
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s.ft('t2')
    s['power',:zero_crossing+1] *= -1
    #}}}
    #{{{plotting enhancement curve at lowest and highest powers 
    idx_maxpower = np.argmax(s.getaxis('power'))
    if fl is not None:
        fl.next('compare highest power to no power')
        fl.plot(s['power',0])
        fl.plot(s['power',idx_maxpower])
    #}}}
    #{{{plotting full enhancement curve
        fl.next('full enhancement curve')
        fl.plot(s)
    s.set_units('power','mW')
    if fl is not None:
        fl.next('real(E(p))')
        fl.image(as_scan_nbr(s.real))
    #}}}
    #{{{plotting enhancement vs power
    enhancement = s['t2':freq_range].sum('t2').real
    enhancement /= enhancement['power',0]
    if fl is not None:
        fl.next('E(p)')
        enhancement.set_units('power','W')
        plt.figure(figsize=(4,4))
        fl.plot((enhancement['power',:idx_maxpower+1]),'ko', human_units=False)
        fl.plot((enhancement['power',idx_maxpower+1:]),'ro', human_units=False)
        plt.title('%s'%searchstr)
        plt.ylabel('Enhancement')
    show()
    return enhancement,idx_maxpower
