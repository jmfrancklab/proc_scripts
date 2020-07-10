from pyspecdata import *
from scipy.optimize import basinhopping
from proc_scripts import *
from proc_scripts import postproc_dict
fl = figlist_var()
logger = init_logging('info')
for searchstr, exp_type, nodename, postproc, label_str in [
        #('200221_CPMG_TEMPOLgel_3p0_1','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        ('200221_CPMG_TEMPOLgel_2p9_1','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200304_CPMG_2p6_1','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_3p5_2','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_3p6_2','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_3p7_2','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_3p7_3','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_3p8_2','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_3p9_2','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        #('200305_CPMG_4p0_1','test_equip','signal','spincore_CPMG_v1','deadtime=5'),
        ]:
    s =  find_file(searchstr, exp_type=exp_type,
            expno=nodename, postproc=postproc, lookup=postproc_dict)
    #{{{ centering CPMG echo
    s = center_CPMG_echo(s)
    fl.next('centered echo')
    fl.image(s)
    #{{{select echo decay fit function
    s.ift('t2')
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    data = s['t2':(0,None)].sum('t2')
    fl.next('Echo decay')
    fl.plot(data,'o')
    print("starting T2 curve")
    f = fitdata(data.real)
    M0,R2,tE = sympy.symbols("M_0 R_2 tE", real=True)
    f.functional_form = M0*sympy.exp(-tE*R2)
    fl.next('T2 test')
    fl.plot(f,'o',label=f.name())
    f.fit()
    fl.plot(f.eval(100),label='%s fit'%f.name())
    text(0.75,0.25, f.latex(),transform=gca().transAxes, size='large',
            horizontalalignment='center', color= 'k')
    print("output",f.output())
    print("latex",f.latex())
    T2 = 1./f.output('R_2')
    fl.show();quit()
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
