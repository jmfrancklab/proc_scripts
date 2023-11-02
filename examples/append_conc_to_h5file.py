from pyspecdata import *
from pyspecProcScripts import *
import h5py
import os, time, h5py
filename = '220126_Ras_M67R1a_capProbe'
actual_conc = 72e-6 #M
h5file = find_file(f"{filename}.h5", exp_type = 'ODNP_NMR_comp/ODNP', expno = 'enhancement')
thisdict = h5file.get_prop('acq_params')
thisdict['concentration'] = actual_conc
h5file.set_prop('acq_params',thisdict)
h5file.name('enhancement')
h5file.hdf5_write(filename + "_rewrite.h5", directory = getDATADIR(exp_type="ODNP_NMR_comp/ODNP"))
quit()

