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

signal_pathway = {'ph1':1,'ph2':-2}
# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
for searchstr,exp_type,nodename,postproc,freq_range,t_range in [
        ["210317_TEMPOL10mM_DNP_cap_probe_1", 'ODNP_NMR_comp', 'signal',
            'spincore_ODNP_v1', (-10000,9000),(None,0.065)]
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
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    #fl.show();quit()
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
    s.ift('t2') # inverse fourier transform into time domain
    logger.debug(strm("THIS IS THE SHAPE"))
    logger.debug(strm(ndshape(s)))
    fl.next('time domain')
    fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    #fl.show();quit()
    s = slice_FID_from_echo(s,max_t=t_range[-1],fl=None)    # visualize time domain after filtering and phasing
    fl.next('FID sliced t domain')
    fl.image(as_scan_nbr(s))
    #{{{apodizing and zero fill
    print(s.dimlabels)
    s.ft('t2')
    fl.next('FID sliced freq domain')
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
    #fl.basename = 'correlation subroutine -- after zero crossing'
    #print(ndshape(s))
    #opt_shift,sigma = correl_align(s['power',zero_crossing+1:],indirect_dim='power',
    #        ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
    #        sigma=50)
    #s.ift('t2')
    #s['power',zero_crossing+1:] *= np.exp(-1j*2*pi*opt_shift*s.fromaxis('t2'))
    #s.ft('t2')
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
    fl.next('apodize and zero fill')
    R = 5.0/(t_range[-1]) # assume full decay by end time
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
    fl.next('1.25 mM TEMPOL ODNP')
    enhancement = s['t2':freq_range].sum('t2').real
    enhancement /= enhancement['power',0]
    enhancement.set_units('power','W')
    plt.figure(figsize=(4,4))
    fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
    fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)
    plt.title('1.25 mM TEMPOL')
    plt.ylabel('Enhancement')
    line = enhancement['power',:idx_maxpower+1]
    x = line.fromaxis('power')
    f=fitdata(line)
    A,B,C,power = symbols("A B  C power")
    f.functional_form = A*s_exp(-B*power)+C
    f.fit()
    fl.plot(f.eval(100),label='fit')
    plt.text(0.75,0.25,f.latex(),transform=plt.gca().transAxes,size='large',
            horizontalalignment='center')
    plt.show()
    init_printing(use_unicode=True)
    print(limit((4.92*s_exp(-5.9*power))-4.59,power,8734))
    quit()
    T1p = nddata(r_[1.298857,1.401268,1.478219,1.568843,1.638806,1.692369],[-1],
            ['power']).setaxis('power',r_[0.001,0.502,1.0,1.58,2.0,2.5])
    fl.next(r'$T_{1}$(p) for 1.25 mM TEMPOL')
    fl.plot(T1p,'o')
    R1w = 1/2.5
    R1p = nddata(r_[0.769908,0.713639,0.67649,0.637412,0.6102,0.59088],[-1],
            ['power']).setaxis('power',r_[0.001,0.502,1.0,1.5,2.0,2.5])
    #{{{making Flinear and fitting
    Flinear = (R1p - R1p['power':0] + R1w)
    polyorder = 3
    coeff,_ = Flinear.polyfit('power',order=polyorder)
    power = nddata(np.linspace(0,R1p.getaxis('power')[-1],25),'power')
    Flinear_fine = 0
    for j in range(polyorder + 1):
        Flinear_fine += coeff[j] * power **j
    fl.next('Flinear')
    fl.plot(Flinear,'o',label='Flinear')
    fl.plot(Flinear_fine,label='Flinear_fine')
    plt.title('polynomial fit of linear equation')
    plt.ylabel("$F_{linear}$")
    fl.next('R1p vs power')
    R1p_fine = (Flinear_fine) + R1p['power':0]-R1w
    lit = nddata(r_[0.85,0.81,0.77,0.745,0.725,0.69],[-1],
            ['power']).setaxis('power',r_[0.001,0.502,1.0,1.5,2.0,2.5])
    fl.plot(lit,'o',label='literature values')
    fl.plot(R1p,"x")
    fl.plot(R1p_fine)
    plt.title("relaxation rates")
    plt.ylabel("$R_1(p)$")
    #{{{plotting without correcting for heating
    ksigs_noT = (1-(enhancement['power',:idx_maxpower+1]))/(659.33*0.00125*T1p['power':0])
    ksigs_noT_max = (1-(-42))/(659.33*0.00125)
    fl.next(r'ksigs_noT vs power')
    fl.plot(ksigs_noT,'o',label='NOT corrected for heating')
    #}}}
    T1p_fine = R1p_fine**-1
    #{{{plotting with correction for heating
    ksigs_T=(1-(enhancement['power',:idx_maxpower+1]))/(659.33*0.00125*T1p_fine)
    fl.plot(ksigs_T,'--',label='with heating correction')
    plt.title('ksigmas(p) vs Power')
    plt.ylabel('ksigmas(p)')
    fl.show();quit()

        #}}}

