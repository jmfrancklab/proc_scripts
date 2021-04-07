from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
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

# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
for searchstr,exp_type,nodename,postproc,freq_range,t_range in [
        ["210325_TEMPOL_10mM_cap_probe_DNP_1", 'ODNP_NMR_comp', 'signal',
            'spincore_ODNP_v1', (-500,100),(None,0.03)]
        #["201203_4AT10mM_DNP_cap_probe_1",'ODNP_NMR_comp','signal',
        #    'spincore_ODNP_v1', (-5000,5000),0.06]
        ]:
    fl.basename = searchstr
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict,fl=fl)
    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range) # visualize the frequency limits
    rcParams.update({
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.9),
        "savefig.facecolor": (1.0,1.0,1.0,0.0),
        })
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
    s.ift('t2') # inverse fourier transform into time domain
    logger.debug(strm("THIS IS THE SHAPE"))
    logger.debug(strm(ndshape(s)))
    fl.next('time domain')
    fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    s = slice_FID_from_echo(s,max_t=t_range[-1],fl=None)    # visualize time domain after filtering and phasing
    #{{{apodizing and zero fill
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
    fl.next('10 mM TEMPOL ODNP Enhancement Curve')
    enhancement = s['t2':freq_range].sum('t2').real
    enhancement /= enhancement['power',0]
    enhancement.set_units('power','W')
    plt.figure(figsize=(4,4))
    fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
    fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)
    plt.title('10 mM TEMPOL Enhancement Curve')
    plt.ylabel('Enhancement')
    plt.show()
    T1p = nddata(r_[0.31,0.28,0.29,0.30,0.32,0.33,0.35],[-1],
            ['power']).setaxis('power',r_[0,0.25,0.5,1,1.5,2,2.5])
    fl.next(r'$T_{1}$(p) for 10 mM TEMPOL')
    fl.plot(T1p,'o')
    R1w = 1/2.8
    R1p = nddata(r_[0.95,1.00,1.04,1.10,1.16,1.19,1.24],[-1],
            ['power']).setaxis('power',r_[0.001,0.251,0.501,1.0,1.58,2.0,2.5])
    #{{{making Flinear and fitting
    Flinear = (R1p - R1p['power':0] + R1w) ** -1
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
    R1p_fine = (Flinear_fine ** -1) + R1p['power':0]-R1w
    fl.plot(R1p,"x")
    fl.plot(R1p_fine)
    plt.title("relaxation rates")
    plt.ylabel("$R_1(p)$")
    plt.show()
    #{{{plotting without correcting for heating
    ksigs_noT = (1-(enhancement['power',:idx_maxpower+1]))/(659.33*0.0025)
    ksigs_noT_max = (1-(-76.2))/(659.33*0.0025)
    fl.next(r'ksigs_noT vs power')
    fl.plot(ksigs_noT,'o',label='NOT corrected for heating')
    #}}}
    #{{{plotting with correction for heating
    ksigs_T=((1-(enhancement['power',:idx_maxpower+1]))/(659.33*0.0025))*R1p_fine
    fl.plot(ksigs_T,'--',label='with heating correction')
    plt.title('ksigmas(p) vs Power')
    plt.ylabel('ksigmas(p)')
    fl.show();quit()

        #}}}

