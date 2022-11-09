from pyspecdata import *
from pyspecProcScripts import *
import h5py
fl = figlist_var()
data_target = os.path.normpath(getDATADIR('ODNP_NMR_comp/inv_rec'))
init_logging(level='debug')
for (filename, IR_kwargs,Ep_kwargs, which_exp) in [
        ('220818_pRA174_MSL_500uM_IR',
            dict(IR_f_slice = (50,250),
                IR_nodenames = [('IR_1')],
                W = 6+1.024, 
                log=False, 
                IR_postproc = None, fl=fl),
            dict(Ep_drif_max = 200),
            dict(has_Ep=False,has_IR=True))
        ]:
    all_kwargs = {**IR_kwargs,**which_exp}
    T1p = generate_T1_Ep(filename, **(all_kwargs))
fl.show()    


