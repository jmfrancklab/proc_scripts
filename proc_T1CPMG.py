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
    #{{{select coherence 
    s = s['ph1',0]['ph2',-1]
    fl.next('select coherence')
    fl.image(s)
    #}}}
    #{{{chunk t2 axis into echoes
    num_echoes = int(s.get_prop('acq')['L'][25])
    logger.debug(strm("here is the t2 axis before chunking the echoes",
        s.getaxis('t2')[r_[0,-1]]))
    s.chunk('t2',['echoes','t2'],[num_echoes,-1])
    logger.debug(strm("here is the t2 axis after chunking the echoes",
        s.getaxis('t2')[r_[0,-1]]))
    fl.next('t2 chunked', figsize=(5,20))
    fl.image(s)
    #}}}
    #{{{centering echoes
    centers = []
    for j in range(ndshape(s)['indirect']):
        s_slice = s['indirect',j]
        this_center = find_echo_center(s_slice)
        centers.append(this_center)
    logger.info(centers)
    avg_center = sum(centers)/len(centers)
    s.ft('t2')
    logger.info(ndshape(s))
    s = center_echo(s, avg_center)
    logger.info(ndshape(s))
    fl.next('s centered')
    logger.info(ndshape(s))
    fl.image(s)
    logger.info(ndshape(s))
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
    fl.show();quit()

