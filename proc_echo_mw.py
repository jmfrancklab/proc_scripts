from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.correlation_alignment_ODNP import correl_align
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
<<<<<<< HEAD
def as_scan_nbr(s):
        return s.C.setaxis('power','#').set_units('power','scan #')
                           
T1p = nddata(r_[1.831,1.970,2.116,2.257,2.318],[-1],
        ['power']).setaxis('power',r_[0.001,0.5,1.0,1.5,2.0])
R1w = 1/2.172
R1p = nddata(r_[0.546,0.508,0.473,0.443,0.431],[-1],
        ['power']).setaxis('power',r_[0.001,0.5,1.0,1.5,2.0])

=======
T1p = nddata(r_[1.831,1.838,1.842,1.97,2.116,2.184,2.257,2.318],[-1],
        ['power']).setaxis('power',r_[0.001,0.05,0.1,0.5,1,1.26,1.58,2])
R1w = 1/2.207
R1p = nddata(r_[0.546,0.544,0.543,0.508,0.473,0.458,0.443,0.431],[-1],
        ['power']).setaxis('power',r_[0.001,0.05,0.1,0.5,1,1.26,1.58,2])
>>>>>>> 73bd74b04ec1601bbc37903868f47768c7eabb34
signal_pathway = {'ph1':1,'ph2':-2}
# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be

# leave this as a loop, so you can load multiple files
<<<<<<< HEAD
for searchstr,exp_type,nodename,postproc,freq_range,t_range,nPowers in [
        ["210518_F195R1a_pR_DHPC_ODNP", 'odnp', 'signal',
            'spincore_ODNP_v1', (-300,300), (None,50e-3),23]
=======
for searchstr,exp_type,nodename,postproc,freq_range,t_range in [
        ["210517_4OHTempo_TempControl_probe_DNP_1", 'ODNP_NMR_comp', 'signal',
            'spincore_ODNP_v1', (-600,100),(None,0.083)]
>>>>>>> 73bd74b04ec1601bbc37903868f47768c7eabb34
        #["201203_4AT10mM_DNP_cap_probe_1",'ODNP_NMR_comp','signal',
        #    'spincore_ODNP_v1', (-5000,5000),0.06]
        ]:
    fl.basename = searchstr
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict,fl=fl)
    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range) # visualize the frequency limits
   # fl.show();quit()
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
<<<<<<< HEAD
    figure()
=======
    s.ift('t2')
>>>>>>> 73bd74b04ec1601bbc37903868f47768c7eabb34
    fl.next('look at data')
    fl.image(s.C.setaxis('power','#'))
    rx_offset_corr = s['t2':(0.065,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft('t2')
    fl.next('')
    fl.image(s.C.setaxis('power','#'))
    s.ift('t2') # inverse fourier transform into time domain
    logger.debug(strm("THIS IS THE SHAPE"))
    logger.debug(strm(ndshape(s)))
    fl.side_by_side('time domain',s,t_range)
    #fl.show();quit()
<<<<<<< HEAD
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))

=======
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway),down_from_max=0.9)
>>>>>>> 73bd74b04ec1601bbc37903868f47768c7eabb34
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
    #R = 5.0/(t_range[-1]) # assume full decay by end time
    #s *= np.exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=2048)
    fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    #}}}
    #{{{select coherence channel in time domain
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s.ft('t2')
    s['power',:zero_crossing+1] *= -1
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

#    plt.figure(figsize=(15,15))
#    for i in range(23):
#        plt.plot(s['power',i].C.data,'-',linewidth=3,alpha=0.5,
#                 label = '%0.2f dBm'%s.getaxis('power')[i])
#        plt.xlim(925,1030)
#        plt.legend(bbox_to_anchor=(1.05,1),loc = 'upper left')
    

    print(s.get_units('power'))
    s.set_units('power','mW')
    fl.next('real(E(p))')
    fl.image(as_scan_nbr(s.real))
    #}}}
    #{{{plotting enhancement vs power
<<<<<<< HEAD
    fl.next('E(p)')
    enhancement = s['t2':freq_range].sum('t2').real
=======
    fl.next('150 uM TEMPOL ODNP')
    #enhancement = s['t2':freq_range].sum('t2').real
    enhancement = s['t2':(-50,50)].sum('t2').real
>>>>>>> 73bd74b04ec1601bbc37903868f47768c7eabb34
    enhancement /= enhancement['power',0]
    enhancement.set_units('power','W')
    #enhancement *= -1
    #enhancement = 1-(enhancement/enhancement.data.max())
    plt.figure(figsize=(4,4))
<<<<<<< HEAD
    fl.plot((1-enhancement['power',:idx_maxpower+1]),'ko', human_units=False)
    fl.plot((1-enhancement['power',idx_maxpower+1:]),'ro', human_units=False)
    plt.ylabel('Enhancement')
    show()
# In[]
=======
    fl.plot((enhancement['power',:idx_maxpower+1]),'ko', human_units=False)
    fl.plot((enhancement['power',idx_maxpower+1:]),'ro', human_units=False)
    fl.show();quit()
    plt.title('1.07 mM TEMPOL')
    plt.ylabel('Enhancement')
    show()
    quit()
>>>>>>> 73bd74b04ec1601bbc37903868f47768c7eabb34
    fl.next(r'$T_{1}$(p) vs power')
    fl.plot(T1p,'o')
    #{{{making Flinear and fitting
    Flinear = ((R1p - R1p['power':0.001] + R1w)**-1)
    print(Flinear)
    polyorder = 3
    coeff,_ = Flinear.polyfit('power',order=polyorder)
    power = nddata(np.linspace(0,R1p.getaxis('power')[-1],nPowers),'power')
    #power = enhancement['power',:idx_maxpower+1].fromaxis('power')    
    Flinear_fine = 0
    for j in range(polyorder + 1):
        Flinear_fine += coeff[j] * power **j
    fl.next('Flinear')
    Flinear.set_units('power',enhancement.get_units('power'))
    fl.plot(Flinear,'o',label='Flinear')
    Flinear_fine.set_units('power',enhancement.get_units('power'))
    fl.plot(Flinear_fine,label='Flinear_fine')
    plt.title('polynomial fit of linear equation')
    plt.ylabel("$F_{linear}$")
    #fl.show();quit()
    fl.next('R1p vs power')
    R1p_fine = ((Flinear_fine)**-1) + R1p['power':0.001]-R1w
    print("UNITS OF R1P:",R1p.get_units('power'))
    fl.plot(R1p,"x")
    R1p_fine.set_units('power',R1p.get_units('power'))
    print("UNITS OF R1P_FINE:",R1p_fine.get_units('power'))
    fl.plot(R1p_fine)
    plt.title("relaxation rates")
    plt.ylabel("$R_1(p)$")
    #{{{plotting without correcting for heating
    ksigs_T=(0.0015167/0.000427)*(1-enhancement['power',:idx_maxpower+1])*(R1p_fine)
    #ksigs_noT = (0.0015167/0.006)*((enhancement['power',:idx_maxpower+1])*(T1p['power':0]**-1))
    fl.next('ksig_smax for 500 uM TEMPOL')
    ksigs_T.set_units('power','mW')
    print("UNITS OF KSIGS_t:",ksigs_T.get_units('power'))
    #ksigs_noT.setaxis('power',ksigs_T.getaxis('power')),
    #ksigs_noT.set_units('power','mW')
    #fl.plot(ksigs_noT,color='r',label='NOT corrected for heating')
    #}}}
    #{{{plotting with correction for heating
    x = enhancement['power',:idx_maxpower+1].fromaxis('power')
    fitting_line = fitdata(ksigs_T)#['power':(0.05,None)])
    k,p_half,power = symbols("k, p_half, power",real=True)
    fitting_line.functional_form = (k*power)/(p_half+power)
    fitting_line.fit()
    fl.plot(ksigs_T,'o',label='with heating correction',human_units=False)
    fl.plot(ksigs_T.imag,'o',label='imaginary',human_units=False)
    print("UNITS OF KSIGS_t IMAG:",ksigs_T.imag.get_units('power'))
    #fl.plot(ksigs_T.imag,'o',label='imaginary')
    fit = fitting_line.eval(25)
    print("fit units:",fit.get_units('power'))
    fit.set_units('power','mW')
    print("fit units after changing:", fit.get_units('power'))
    fl.plot(fit,label='fit',human_units=False)
    plt.text(0.75, 0.25, fitting_line.latex(), transform=plt.gca().transAxes,size='large',
            horizontalalignment='center',color='k')
    #ax = plt.gca()
    plt.title('ksigmas(p) vs Power')
    plt.ylabel('ksigmas(p)')
    fl.show()

        #}}}

# In[]