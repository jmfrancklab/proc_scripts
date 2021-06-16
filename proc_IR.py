from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from pyspecProcScripts.third_level.process_IR import process_IR
fl = fl_mod()
# {{{ input parameters
save_npz = False
#}}}
coherence_pathway = {'ph1':0,'ph2':1}
for thisfile,exp_type,nodename,postproc,f_range,t_range,clock_correction,IR,ILT in [
       ('210615_S175R1a_pR_DDM_ODNP','odnp',
           'FIR_0dBm','spincore_IR_v1',
           (-0.5e3,0.5e3),(None,50e-3),True,False,False),
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    myslice = s['t2':f_range]
    mysgn = determine_sign(select_pathway(myslice, coherence_pathway), fl=fl)
    T1 = process_IR(s,label=thisfile,W=7,f_range=f_range,t_range=t_range,
            clock_correction=clock_correction,IR=IR,flip=True,sign=mysgn,fl=fl) 
    fl.show()
