"""
Edit the concentration parameter of a saved HDF5 file
=====================================================

Loads in a previously saved HDF5 dataset. If the attribute
'concentration' does not exist in the acq_params, this will
add the concentration defined by 'actual_conc' into acq_params
under 'concentration'. If the attribute already exists it will
simply replace the current concentration with the value of
'actual_conc'.
"""

from pyspecdata import *
import h5py, os

data_info = dict(
    filename="260526_hydroxytempo_ODNP_2.h5",  # file that is being edited
    file_location="B27/ODNP",  # location of file
    actual_conc=27e-3,
)  # concentration of the dataset in M
h5 = search_filename(
    data_info["filename"], exp_type=data_info["file_location"], unique=True
)
with h5py.File(os.path.normpath(h5), "r+") as thisfile:
    for node in thisfile:
        if "other_info/acq_params" not in thisfile[node]:
            continue
        acq_params = thisfile[node]["other_info"]["acq_params"]
        acq_params.attrs["concentration"] = data_info["actual_conc"]
    thisfile.close()
