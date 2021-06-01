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
def process_enhancement(s, searchstr='', signal_pathway = {'ph1':1,'ph2':0},
        excluded_pathways = [(0,0),(0,1)], freq_range=(None,None),
        t_range=(0,0.083),fl=None):
    #if fl is not None:
    #    fl.basename = searchstr
    ##    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
    #            s,freq_range) # visualize the frequency limits
    s.ift('t2')
    s.reorder(['ph1','ph2','power','t2'])
    #if fl is not None:
    #    fl.basename = searchstr
    #    fl.side_by_side('show time limits',s,t_range)
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s.ift(['ph1','ph2'])   
    t_start = t_range[-1]/4
    t_start *= 3
    rx_offset_corr = s['t2':(t_start,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    s.ft(['ph1','ph2'])
    zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('power',raw_index=True).item()
    s = s['t2':freq_range]    
    s.ift('t2') # inverse fourier transform into time domain
    best_shift,max_shift = hermitian_function_test(select_pathway(s,signal_pathway).C.convolve('t2',0.01))
    best_shift = 0.33e-3
    s.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
    s.ft('t2')
    logger.info(strm("applying zeroth order correction"))
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',1]['ph2',0] = 1
    phasing.ift(['ph1','ph2'])
    s /= phasing
    ph0 = s.C.sum('t2')
    ph0 /= abs(ph0)
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
    s.setaxis('power',power_axis_W)
    s.ift(['ph1','ph2'])
    #if fl is not None:
    #    fl.next('phased-frequency domain')
    #    fl.image(as_scan_nbr(s))
    opt_shift,sigma = correl_align(s,indirect_dim='power',
            ph1_selection=1,ph2_selection=0,sigma=50)
    s.ift('t2')
    s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename= None
    if fl is not None:
        fl.next(r'after correlation, $\varphi$ domain')
        fl.image(as_scan_nbr(s))
    #fl.show();quit()
    s *= phasing
    s.ift('t2')
    s.ft(['ph1','ph2'])
    #if fl is not None:
    #    fl.next('after correlation alignment FTed ph')
    #    fl.image(as_scan_nbr(s))
    s.reorder(['ph1','ph2','power','t2'])
    #if fl is not None:
    #    fl.next('after correlation -- time domain')
    #    fl.image(as_scan_nbr(s))
    s.ft('t2')    
    #if fl is not None:
    #    fl.next('after correlation -- frequency domain')
    #    fl.image(as_scan_nbr(s))
    s.ift('t2')
    d=s.C
    
    #d['power',zero_crossing] *= -1
    #d *= -1
    d.ft('t2')
    d.ift('t2')
    d = d['t2':(0,None)]
    d['t2':0] *= 0.5
    d.ft('t2')
    # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
    error_pathway = (set(((j,k) for j in range(ndshape(d)['ph1']) for k in range(ndshape(d)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    # }}}
    d_,frq_slice = integral_w_errors(d,signal_pathway,error_pathway,
            indirect='power', fl=fl, return_frq_slice=True)
    #if fl is not None:
    #    fl.next('dummy')
    #    with figlist_var() as fl:
    #        fl.next('extra')
    #        left_pad,bottom_pad,width_pad,top_pad = DCCT(s,fl.next('dummy'),just_2D=True)
    #        pass_frq_slice=True
    #        frq_slice=frq_slice
    #        fl.next('diagnostic 1D plot')
    #        fl.plot(d['power',:]['ph1',signal_pathway['ph1']]['ph2',signal_pathway['ph2']].real,alpha=0.4)
    #        axvline(x=frq_slice[0],c='k',linestyle=':',alpha=0.8)
    #        axvline(x=frq_slice[-1],c='k',linestyle=':',alpha=0.8)
    #        fl.next('')
    #        x = d.getaxis('t2')
    #        dx = x[1]-x[0]
    #        y = d.getaxis('power')
    #        dy = y[1]-y[0]
    #        start_y = d.getaxis('power')[0]
    #        stop_y = d.getaxis('power')[-1]+25
    #        figure()
    #        ax = plt.axes([left_pad,bottom_pad,width_pad,top_pad])
    #        yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,10])
    #        majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=2, steps=[1,10])
    #        minorLocator = lambda: mticker.AutoMinorLocator(n=4)
    #        ax.xaxis.set_major_locator(majorLocator())
    #        ax.xaxis.set_minor_locator(minorLocator())
    #        ax.set_ylabel(None)
    #        fl.image(d['ph1',signal_pathway['ph1']]['ph2',signal_pathway['ph2']].real.run(complex128).C.setaxis(
#'power','#').set_units('power','scan #'),black=False)
    #        x1 = x[0]-dx
    #        y1 = start_y-dy
    #        x2 = frq_slice[-1]+dx
    #        tall = 500#(y[-1]+dy)-(y[0]-dy) 
    #        wide = frq_slice[0]-x1
    #        wide2 = (x[-1]+dx)-x2
    #        p = patches.Rectangle((x1,y1),wide,tall,angle=0.0,linewidth=1,fill=None,
    #                hatch='//',ec='k')
    #        ax.add_patch(p)
    #        q = patches.Rectangle((x2,y1),wide2,tall,angle=0.0,linewidth=1,
    #                fill=None,hatch='//',ec='k')
    #        ax.add_patch(q)
    d = d_.C
    idx_maxpower = np.argmax(s.getaxis('power'))
    if fl is not None:
        fl.next('full enhancement curve')
        fl.plot(d)
    d /= max(d.data)
    d.setaxis('power',power_axis_W)
    thiscolor = next(thesecolors)
    if fl is not None:
        fl.next('E(p)')
        fl.plot(d['power',:-3], 'ko', capsize=6, alpha=0.3)
        fl.plot(d['power',-3:],'ro',capsize=6, alpha=0.3)
    enhancement = d
    show()
    return enhancement,idx_maxpower
