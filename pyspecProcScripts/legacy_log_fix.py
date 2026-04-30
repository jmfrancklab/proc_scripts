import datetime
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pyspecdata as psd

import pyspecProcScripts as prscr


def attach_log_data_from_file(
    s,
    filename,
    exp_type,
    node_name="log",
    hdf_repair=None,
):
    """Attach legacy HDF log data as an nddata property.

    Because the format of the log has not changed **only** its location
    inside the HDF5, this is the only code that we need to maintain
    backwards-compatibility.
    """
    filename = psd.search_filename(
        re.escape(filename), exp_type=exp_type, unique=True
    )
    with h5py.File(filename, "r") as f:
        if hdf_repair is None:
            thislog = prscr.logobj.from_group(f[node_name])
        else:
            thislog = prscr.logobj.from_group(hdf_repair(f[node_name]))
        # TODO ☐: previously there was extra code here, which could be
        #         unnecessary.  I'm not sure what the reason for this
        #         is.  Looking at test_logobj.py, which is in FLInst,
        #         from_group should give exactly what we need
    s.set_prop("log", thislog)
    return s
