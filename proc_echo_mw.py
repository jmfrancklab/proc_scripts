from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
rcParams["savefig.transparent"] = True
logger = init_logging("info")
fl = fl_mod()
t2 = symbols('t2')

# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
# leave this as a loop, so you can load multiple files
for searchstr,exp_type,nodename,postproc,freq_range,max_t in [
        ["201208_Ni_sol_probe_6", 'ODNP_NMR_comp', 'signal',
            'spincore_ODNP_v1', (-300,300), 0.055]
        ]:
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict)
    quit()
    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range) # visualize the frequency limits
    s = s['t2':freq_range] # slice out the frequency range along t2 axis
    s.ift('t2') # inverse fourier transform into time domain
    logger.debug(strm("THIS IS THE SHAPE"))
    logger.debug(strm(ndshape(s)))
    s = slice_FID_from_echo(s,max_t=max_t,fl=fl)    # visualize time domain after filtering and phasing
    #{{{apodizing and zero fill
    fl.next('apodize and zero fill')
    R = 5.0/(max_t) # assume full decay by end time
    s *= exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1024)
    fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    print(s.getaxis('power'))
    #}}}
    #{{{select coherence channel in time domain
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s.ft('t2')
    
    #}}}
    #{{{plotting enhancement curve at lowest and highest powers 
    fl.next('compare highest power to no power')
    idx_maxpower = argmax(s.getaxis('power'))
    fl.plot(s['power',0])
    fl.plot(s['power',idx_maxpower])
    #}}}
    #{{{plotting full enhancement curve
    fl.next('full enhancement curve')
    fl.plot(s)
    #fl.show();quit()
    #}}}
    #{{{plotting enhancement vs power
    fl.next('enhancement for 695 uM 4-AT')
    enhancement = s['t2':(-400,400)].sum('t2').real
    enhancement /= enhancement['power',0]
    enhancement.set_units('power','W')
    fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
    fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)
    ylabel('Enhancement')
    print(ndshape(enhancement))
    enhancement = enhancement['power':4]
    print(enhancement) 
    #}}}
fl.show();quit()
