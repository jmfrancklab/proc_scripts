from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
fl = fl_mod()
for searchstr,exp_type,postproc in [
        ["201118_1mM4AT",'ESR','ESR_linewidth']
        ]:
    s = find_file(searchstr+'.DSC',exp_type=exp_type,postproc=postproc,
            lookup=postproc_dict)
    fl.next('linewidth for 100 uM 4AT')
    fl.plot(s)
    
    fl.show();quit()

