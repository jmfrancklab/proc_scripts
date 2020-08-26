from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import fl_mod,find_echo_center, center_echo
from proc_scripts import postproc_dict
from sympy import symbols
from proc_scripts.fitting import decay
logger = init_logging("debug")
logger.info("this is a test")
logger.info("this is a test debug")

fl = fl_mod()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
clock_correction = 0
filter_bandwidth = 5e3
t2 = symbols('t2')
# }}}
for searchstr, exp_type, nodename in [
        ('w8_200731','NMR_Data_AG',5)
        #('200303','T1CPMG_AER')
        ]:
    s = find_file(searchstr,exp_type=exp_type,
            expno=nodename,lookup=postproc_dict, fl=fl)
    #{{{centering echoes
    s.ft('t2',shift=True)
    s.ift('t2')
    centers = []
    for j in range(ndshape(s)['indirect']):
        s_slice = s['indirect',j]
        this_center = find_echo_center(s_slice,fl=fl)
        centers.append(this_center)
    logger.info(centers)
    avg_center = sum(centers)/len(centers)
    logger.info(strm("averaged center is",avg_center))
    print("here is the t2 axis before centering",
        s.getaxis('t2')[r_[0,-1]])
    s = center_echo(s, avg_center)
    fl.next('s centered')
    fl.image(s)
    fl.show();quit()
    #}}}
    #}}}
    #{{{fitting decay function
    s.ft('t2')
    fl.next('freq domain')
    fl.image(s)
    fl.show();quit()
    f,T2 = decay(s, (None,None), indirect='echoes')
    fl.plot_curve(f,'T2 relaxation decay')
    how();quit()

