from pyspecdata import *
from scipy.optimize import basinhopping
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.fitting import decay
fl = fl_mod()
logger = init_logging('debug')
dwdel1=6.5e-6
TD=8
tau_extra=20e-6
for searchstr, exp_type, nodename, postproc, label_str, f_range, spincore in [
        #('w8_200917','test_equip',6,'ag_CPMG_strob','water loading 8',(-500,500),False),
        #('freeSL_201001','test_equip',7,'ag_CPMG_strob','free SL',(-500,500),False),
        #('200221_CPMG_TEMPOLgel_2p9_1','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500),True),
        #('w8_200731','test_equip',3,'ag_CPMG_strob','water loading 8',(-500,500),False),
        ('freeSL_201007','test_equip',3,'ag_CPMG_strob','free SL',(-500,500),False),
        #('200305_CPMG_3p5_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p6_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p7_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p7_3','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p8_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_3p9_2','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        #('200305_CPMG_4p0_1','test_equip','signal','spincore_CPMG_v1','deadtime=5',(-500,500)),
        ]:
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename, postproc=postproc, lookup=postproc_dict, fl=fl)
    fl.show();quit()
    s.ift('t2')
    if spincore:
        s.reorder('nScans',first=True)
        s = s['ph1',1]
        s.mean('nScans')
        s.reorder('t2',first=True)
    #{{{ centering CPMG echo
    center = find_echo_center(s)
    s = center_echo(s,center,fl=fl)
    logger.debug(strm(ndshape(s)))
    fl.next('centered echo')
    fl.image(s)
    #{{{select echo decay fit function
    s.ft('t2')
    fl.next('before fitting')
    fl.image(s)
    #fl.show();quit()
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
