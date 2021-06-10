from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.third_level.process_IR import process_IR
fl = fl_mod()
# {{{ input parameters
save_npz = False
#}}}
coherence_pathway = {'ph1':0,'ph2':1}
for thisfile,exp_type,nodename,postproc,f_range,t_range,clock_correction,IR,ILT in [
       ('210610_water_TempCont_probe_DNP_1','ODNP_NMR_comp/ODNP',
           'FIR_34dBm','spincore_IR_v1',
           (0.1e3,0.345e3),(None,50e-3),True,False,False),
        # ('210322_water_control_FIR_noPower','inv_rec','signal','spincore_IR_v1',
        #    (-0.09e3,-0.06e3),(None,83e-3),False,False),
        #('210517_4OHTempo_TempControl_probe_FIR_34dBm.','odnp_nmr_comp/inv_rec','signal','spincore_IR_v1',
        #    (-0.4e3,0.4e3),(None,20e-3),False,False),
        #('w3_201111','test_equip',2,'ag_IR2H',(-600,600),(0,None),True)
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    myslice = s['t2':f_range]
    mysgn = determine_sign(select_pathway(myslice, coherence_pathway), fl=fl)
    T1 = process_IR(s,label=thisfile,W=7,f_range=f_range,t_range=t_range,
            clock_correction=clock_correction,IR=IR,flip=True,sign=mysgn,fl=fl) 
    fl.show()
