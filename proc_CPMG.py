from pyspecdata import *
from scipy.optimize import basinhopping
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.fitting import decay
fl = fl_mod()
logger = init_logging('info')
for searchstr, exp_type, nodename, postproc, label_str, f_range in [
        #('200221_CPMG_TEMPOLgel_2p9_1','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200304_CPMG_2p6_1','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p5_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p6_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        ('200305_CPMG_3p7_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p7_3','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p8_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p9_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_4p0_1','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        ]:
    s =  find_file(searchstr, exp_type=exp_type,
            expno=nodename, postproc=postproc, lookup=postproc_dict)
    #{{{ centering CPMG echo
    s = center_CPMG_echo(s)
    fl.next('centered echo')
    fl.image(s)
    #{{{select echo decay fit function
    s.ft('t2')
    f,T2 = decay(s, f_range, indirect='tE')
    fl.plot_curve(f,'T2 relaxation decay')
    #}}}
    #{{{saving figure
    save_fig = False
    if save_fig:
        savefig('20200108_CPMG_trials.png',
                transparent=True,
                bbox_inches='tight',
                pad_inches=0,
                legend=True)
    fl.show()
    #}}}
