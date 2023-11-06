from pyspecdata import *
import h5py
filename = '220126_Ras_M67R1a_capProbe.h5'
actual_conc = 72e-6 #M
h5 = search_filename(f"{filename}" , exp_type = 'ODNP_NMR_comp/ODNP')
h5 = ['/'.join(path.split('\\')) for path in h5]
with h5py.File(h5[0],'r+') as thisfile:
    acq_params = thisfile['enhancement']['other_info']['acq_params']
    acq_params.attrs['concentration'] = actual_conc
    thisfile.close()

