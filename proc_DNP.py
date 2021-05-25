from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.third_level.process_IR import process_IR
fl = fl_mod()
# {{{ input parameters
save_npz = False
#}}}
T1_list = []
for thisfile,exp_type,nodename,postproc,f_range,t_range,IR,ILT in [
       ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','0p001_Watts','spincore_IR_v1',
           (-0.5e3,1e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','0p5_Watts','spincore_IR_v1',
           (-0.5e3,1e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','1_Watts','spincore_IR_v1',
           (-0.5e3,1e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','1p5_Watts','spincore_IR_v1',
           (-0.5e3,1e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','2_Watts','spincore_IR_v1',
           (-0.5e3,1e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','2p5_Watts','spincore_IR_v1',
           (-0.5e3,1e3),(None,83e-3),True,False),

        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    #fl.show();quit()
    T1 = process_IR(s,label=thisfile,W=5,f_range=f_range,fl=fl) 
    T1_list.append(T1)
print(T1_list)
quit()

