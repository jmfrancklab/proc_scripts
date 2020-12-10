from pyspecdata import *
from proc_scripts import postproc_dict, fl_mod
from numpy import *
fl = fl_mod()
for searchstr,exp_type,nodename,postproc in [
    ['201209_Ni_sol_probe_var_tau','var_tau','var_tau','spincore_var_tau_v1']
    ]:
    s = find_file(searchstr,exp_type=exp_type,expno=nodename,postproc=postproc,
            lookup=postproc_dict,fl=fl)
    fl.show();quit()

