from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols, Symbol, latex,limit,init_printing
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
fl = fl_mod()
t2 = symbols('t2')
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
T1p = nddata(r_[0.436,0.465,0.482,0.507,0.527],[-1],
        ['power']).setaxis('power',r_[0.001,0.5,1.0,1.5,2.0])
R1w = 1/2.172
R1p = nddata(r_[2.292,2.150,2.073,1.971,1.897],[-1],
        ['power']).setaxis('power',r_[0.001,0.5,1.0,1.5,2.0])
signal_pathway = {'ph1':1,'ph2':-2}
# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
for searchstr,exp_type,nodename,postproc,freq_range,t_range in [
        ["210414_TEMPOL6mM_DNP_cap_probe_1", 'ODNP_NMR_comp', 'signal',
            'spincore_ODNP_v1', (-11000,11000),(None,0.083)]
        #["201203_4AT10mM_DNP_cap_probe_1",'ODNP_NMR_comp','signal',
        #    'spincore_ODNP_v1', (-5000,5000),0.06]
        ]:
    fl.basename = searchstr
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict,fl=fl)
    def as_scan_nbr(s):
        return s.C.setaxis('power','#').set_units('power','scan #')
    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range) # visualize the frequency limits
    #fl.show();quit()
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
    fl.next('look at data')
    fl.image(s.C.setaxis('power','#'))
    s.ift('t2')
    rx_offset_corr = s['t2':(0.02,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    fl.next('look at data T')
    fl.image(s.C.setaxis('power','#'))
    s.ift('t2') # inverse fourier transform into time domain
    logger.debug(strm("THIS IS THE SHAPE"))
    logger.debug(strm(ndshape(s)))
    fl.side_by_side('time domain',s,t_range)
    #fl.show();quit()
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    coh_slice = select_pathway(s['t2':0],signal_pathway)
    fl.next('before zeroth order phasing correction')
    s.ft('t2')
    fl.image(as_scan_nbr(s))
    s.ift('t2')
    print("applying zeroth order correction")
    s.ift(['ph1','ph2'])
    phasing = s['t2',0].C
    phasing.data *= 0
    phasing.ft(['ph1','ph2'])
    phasing['ph1',1]['ph2',0] = 1
    phasing.ift(['ph1','ph2'])
    fl.next('zeroth order corrected')
    ph0 = s['t2':0]/phasing
    ph0 /= abs(ph0)
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
    s /= ph0
    s.ft(['ph1','ph2'])
    #{{{apodizing and zero fill
    print(s.dimlabels)
    s.ft('t2')
    fl.next('phase corrected')
    fl.image(as_scan_nbr(s))
    #fl.show();quit()
    s.reorder(['ph1','ph2','power','t2'])
    zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin('power',raw_index=True).item()
    logger.info(strm("zero corssing at",zero_crossing))
    s.ift(['ph1','ph2'])
    frq_max = abs(s).argmax('t2')
    s.ift('t2')
    s *= np.exp(-1j*2*pi*frq_max*s.fromaxis('t2'))
    s.ft('t2')
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
    print(ndshape(s))
    opt_shift,sigma = correl_align(s['power',zero_crossing+1:],indirect_dim='power',
            ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
            sigma=50)
    s.ift('t2')
    s['power',zero_crossing+1:] *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    s.ft('t2')
    fl.basename = None
    fl.next(r'after correlation, $\varphi$ domain')
    fl.image(as_scan_nbr(s))
    s.ft(['ph1','ph2'])
    fl.next('after correlation alignment FTed ph')
    fl.image(as_scan_nbr(s))
    s.reorder(['ph1','ph2','power','t2'])
    fl.next('after correlation -- frequency domain')
    fl.image(as_scan_nbr(s))
    s.ift('t2')
    fl.next('after correlation -- time domain')
    fl.image(as_scan_nbr(s))
    #fl.show();quit()
    s = s['t2':(0,t_range[-1])]
    s['t2':0] *= 0.5
    if select_pathway(s['t2':0]['power',0],signal_pathway).real < 0:
        s *= -1
    fl.next('apodized and zero fill')
    R = 5.0/(t_range[-1]) # assume full decay by end time
    s *= np.exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1781)
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
    print(s.get_units('power'))
    s.set_units('power','mW')
    #fl.side_by_side('Real of E(p)',as_scan_nbr(s.real),(-500,300),human_units=False)
    fl.next('real(E(p))')
    fl.image(as_scan_nbr(s.real))
    #print(ndshape(s))
    #s_min = s.real.data.min()
    #s_max = s.real.data.max()
    #d = (s.real-s_min)/(s_max-s_min)
    #fl.next('normalized real')
    #fl.image(d.C.setaxis(
#'power','#').set_units('power','scan #'))
    #}}}
    #{{{plotting enhancement vs power
    fl.next('150 uM TEMPOL ODNP')
    enhancement = s['t2':freq_range].sum('t2').real
    enhancement /= enhancement['power',0]
    enhancement.set_units('power','W')
    plt.figure(figsize=(4,4))
    fl.plot((1-enhancement['power',:idx_maxpower+1]),'ko', human_units=False)
    fl.plot((1-enhancement['power',idx_maxpower+1:]),'ro', human_units=False)
    plt.title('150 uM TEMPOL')
    plt.ylabel('Enhancement')
    show()
    fl.next(r'$T_{1}$(p) vs power')
    fl.plot(T1p,'o')
    #{{{making Flinear and fitting
    Flinear = ((R1p - R1p['power':0.001] + R1w)**-1)
    print(Flinear)
    polyorder = 3
    coeff,_ = Flinear.polyfit('power',order=polyorder)
    power = nddata(np.linspace(0,R1p.getaxis('power')[-1],25),'power')
    #power = enhancement['power',:idx_maxpower+1].fromaxis('power')    
    Flinear_fine = 0
    for j in range(polyorder + 1):
        Flinear_fine += coeff[j] * power **j
    fl.next('Flinear')
    fl.plot(Flinear,'o',label='Flinear')
    Flinear_fine.set_units('power',Flinear.get_units('power'))
    fl.plot(Flinear_fine,label='Flinear_fine')
    plt.title('polynomial fit of linear equation')
    plt.ylabel("$F_{linear}$")
    #fl.show();quit()
    fl.next('R1p vs power')
    R1p_fine = ((Flinear_fine)**-1) + R1p['power':0.001]-R1w
    fl.plot(R1p,"x")
    R1p_fine.set_units('power',R1p.get_units('power'))
    fl.plot(R1p_fine)
    plt.title("relaxation rates")
    plt.ylabel("$R_1(p)$")
    #{{{plotting without correcting for heating
    ksigs_T=(0.0015167/0.006)*(enhancement['power',:idx_maxpower+1])*(R1p_fine)
    ksigs_noT = (0.0015167/0.006)*((enhancement['power',:idx_maxpower+1])*(T1p['power':0]**-1))
    fl.next('ksig_smax for 2.25 mM TEMPOL')
    ksigs_noT.setaxis('power',ksigs_T.getaxis('power')),
    ksigs_noT.set_units('power','mW')
    #fl.plot(ksigs_noT,color='r',label='NOT corrected for heating')
    #}}}
    #{{{plotting with correction for heating
    x = ksigs_T.fromaxis('power')
    fitting_line = fitdata(ksigs_T['power':(0.25,None)])
    k,p_half,power = symbols("k, p_half, power",real=True)
    fitting_line.functional_form = (k*power)/(p_half+power)
    fitting_line.fit()
    fl.plot(ksigs_T,'o',label='with heating correction')
    fl.plot(ksigs_T.imag,'o',label='imaginary')
    #fl.plot(ksigs_T.imag,'o',label='imaginary')
    fit = fitting_line.eval(100)
    fit.set_units('power','W')
    fl.plot(fit,label='fit',human_units=False)
    plt.text(0.75, 0.25, fitting_line.latex(), transform=plt.gca().transAxes,size='large',
            horizontalalignment='center',color='k')
    ax = plt.gca()
    plt.title('ksigmas(p) vs Power')
    plt.ylabel('ksigmas(p)')
    fl.show();quit()

        #}}}

