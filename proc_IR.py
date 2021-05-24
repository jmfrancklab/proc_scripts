from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.third_level.process_IR import process_IR
fl = fl_mod()
# {{{ input parameters
save_npz = False
#}}}
for thisfile,exp_type,nodename,postproc,f_range,t_range,IR,ILT in [
       ('210512_water_cap_probe_IR_33dBm','inv_rec','signal','spincore_IR_v1',
           (-0.291e3,0.615e3),(None,83e-3),True,False),
        # ('210322_water_control_FIR_noPower','inv_rec','signal','spincore_IR_v1',
        #    (-0.146e3,-0.06e3),(None,83e-3),False,False),
        #('210517_4OHTempo_TempControl_probe_FIR_34dBm.','odnp_nmr_comp/inv_rec','signal','spincore_IR_v1',
        #    (-0.4e3,0.4e3),(None,20e-3),False,False),
        #('w3_201111','test_equip',2,'ag_IR2H',(-600,600),(0,None),True)
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    process_IR(s,label=thisfile,fl=fl) 

