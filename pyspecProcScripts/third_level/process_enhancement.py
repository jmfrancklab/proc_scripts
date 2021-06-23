from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from pyspecProcScripts.correlation_alignment_ODNP import correl_align
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
def process_enhancement(s, searchstr='', signal_pathway = {'ph1':1},
        excluded_pathways = [(0,0)], freq_range=(None,None),
        t_range=(0,0.083),sign=None,fl=None):
    s *= sign
    #if fl is not None:
    #    print(s.get_units('t2'))
    #    s.set_units('t2','kHz')
    #    fl.push_marker()
    #    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
    #            s,thisrange=freq_range) # visualize the frequency limits
    #fl.show();quit()
    s.ift('t2')
    s.reorder(['ph1','power','t2'])
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s.ift(['ph1'])  
    #{{{ Applying DC offset correction
    t_start = t_range[-1]/4
    t_start *= 3
    rx_offset_corr = s['t2':(t_start,None)]
    rx_offset_corr = rx_offset_corr.data.mean()
    s -= rx_offset_corr
    s.ft('t2')
    s.ft(['ph1'])
    #}}}
    zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('power',raw_index=True).item()
    s = s['t2':freq_range] 
    if fl is not None:
        fl.next('freq_domain before hermitian')
        fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
        s.ift('t2')
        fl.next('time domain before hermitian')
        fl.image(s)
    #{{{Applying phasing corrections
    #s.ift('t2') # inverse fourier transform into time domain
    best_shift,max_shift = hermitian_function_test(select_pathway(s,signal_pathway).C.convolve('t2',0.01))
    best_shift = 0.033e-3
    s.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
    logger.info(strm("applying zeroth order correction"))
    #if fl is not None:
    #    fl.next('After hermitian-t domain')
    #    fl.image(s)
    #    s.ft('t2')
    #    fl.next('after hermitian-freq domain')
    #    fl.image(s)
    #fl.show();quit()    
    s.ift(['ph1'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1'])
    phasing['ph1',1] = 1
    phasing.ift(['ph1'])
    s /= phasing
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    s /= ph0
    s.ft(['ph1'])
    logger.info(strm(s.dimlabels))
    s.ft('t2')
    if fl is not None:
        fl.next('After zeroth order phase correction')
        fl.image(as_scan_nbr(s))
        #s.ift('t2')
        #fl.next('time domain after zeroth order')
        #fl.image(s)
    #fl.show();quit()    
    s.reorder(['ph1','power','t2'])
    logger.info(strm("zero corssing at",zero_crossing))
    #}}}
    #{{{Applying correlation alignment
    s.ift(['ph1'])
    opt_shift,sigma = correl_align(s,indirect_dim='power',
            ph1_selection=1,sigma=0.001)
    s.ift('t2')
    s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename= None
    if fl is not None:
        fl.next(r'after correlation, $\varphi$ domain')
        s.set_units('t2','Hz')
        #s.C.reorder(['power','ph1']).smoosh(['power','ph1'])
        #s /= phasing
        fl.image(s,human_units=False)
    s.ift('t2')
    if fl is not None:
        fl.next('t domain after correlation')
        fl.image(s)
    #fl.show();quit()    
    s.ft(['ph1'])
    if fl is not None:
        fl.next('after correlation alignment FTed ph')
        fl.image(as_scan_nbr(s))
    s.reorder(['ph1','power','t2'])
    if fl is not None:
        fl.next('after correlation -- time domain')
        fl.image(as_scan_nbr(s))
    s.ft('t2')    
    if fl is not None:
        fl.next('after correlation -- frequency domain')
        fl.image(as_scan_nbr(s))
    #}}}
    s.ift('t2')
    d=s.C
    d.ft('t2')
    d.ift('t2')
    d = d['t2':(0,t_range[-1])]
    d['t2':0] *= 0.5
    d.ft('t2')
    if fl is not None:
        fl.next('FID sliced')
        fl.image(d)
        #d.ift('t2')
        #fl.next('FID sliced')
        #fl.image(d,human_units=False)
        #d.ft('t2')
    d *= sign
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    #d.mean('nScans')
    error_pathway = (set(((j) for j in range(ndshape(d)['ph1'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'])]))
    error_pathway = [{'ph1':j} for j in error_pathway]
    # }}}
    #{{{ integrating with error bar calculation
    d_,frq_slice,std = integral_w_errors(d,signal_pathway,error_pathway,
            indirect='power', fl=fl, return_frq_slice=True)
    x = d_.get_error()
    x[:] /= sqrt(2)
    d = d_.C
    #}}}
    
    #{{{Normalizing by max 
    idx_maxpower = np.argmax(s.getaxis('power'))
    d /= max(d.data)
    #}}}
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
    power_axis_W = r_[0,power_axis_W]
    d.setaxis('power',power_axis_W)
    thiscolor = next(thesecolors)
    if fl is not None:
        fl.next('E(p)')
        fl.plot(d['power',:-3], 'ko', capsize=6, alpha=0.3)
        fl.plot(d['power',-3:],'ro',capsize=6, alpha=0.3)
        fl.pop_marker()
    enhancement = d['power',:-3]
    return enhancement,idx_maxpower
