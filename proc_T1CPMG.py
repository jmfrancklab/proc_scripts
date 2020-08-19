from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import fl_mod,center_CPMG_echo
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
        ('w8_200731','test_equip',5)
        #('200303','T1CPMG_AER')
        ]:
    s = find_file(searchstr,exp_type=exp_type,
            expno=nodename,lookup=postproc_dict, fl=fl)
    #{{{select coherence 
    s = s['ph1',0]
    s = s['ph2',-1]
    fl.next('select coherence')
    fl.image(s)
    #}}}
    #{{{chunk t2 axis into echoes
    s.chunk('t2',['echoes','t2'],[128,-1])
    fl.next('t2 chunked', figsize=(5,20))
    fl.image(s)
    print(ndshape(s))
    fl.show();quit() 
    #}}}
    #{{{centering echoes
    s = center_CPMG_echo(s)
    fl.next('centered with center cpmg echo function')
    fl.image(s)
    fl.show();quit()
    #}}}
    #{{{fitting decay function
    s.ft('t2')
    fl.next('freq domain')
    fl.image(s)
    fl.show();quit()
    f,T2 = decay(s, (None,None), indirect='echoes')
    fl.plot_curve(f,'T2 relaxation decay')
    fl.show();quit()

