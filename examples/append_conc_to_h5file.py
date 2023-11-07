"""Edit the concentration parameter of a saved HDF5 file
========================================================

Loads in a previously saved HDF5 dataset. If the attribute
'concentration' does not exist in the acq_params, this will 
add the concentration defined by 'actual_conc' into acq_params
under 'concentration'. If the attribute already exists it will
simply replace the current concentration with the value of 
'actual_conc'.
"""
from pyspecdata import *
import h5py
filename = '220126_Ras_M67R1a_capProbe.h5'
actual_conc = 72e-6 #M
h5 = search_filename(f"{filename}" , exp_type = 'ODNP_NMR_comp/ODNP',unique=True)
print(h5)
h5 = ['/'.join(path.split('\\')) for path in h5]
print(h5)
with h5py.File(h5[0],'r+') as thisfile:
    acq_params = thisfile['enhancement']['other_info']['acq_params']
    acq_params.attrs['concentration'] = actual_conc
    thisfile.close()

