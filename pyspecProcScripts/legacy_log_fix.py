import re
import h5py
import pyspecdata as psd
import pyspecProcScripts as prscr


def attach_log_data_from_file(
    s,
    filename,
    exp_type,
    node_name="log",
):
    """Attach legacy HDF log data as an nddata property.

    Because the format of the log has not changed **only** its location
    inside the HDF5, this is the only code that we need to maintain
    backwards-compatibility.
    This code also solves broken hdf files by reading the log data into
    memory and then writing it back.
    """
    filename = psd.search_filename(
        re.escape(filename), exp_type=exp_type, unique=True
    )
    with h5py.File(filename, "r") as f:
        thislog = prscr.logobj.from_group(f[node_name])
        thislog.total_log = thislog.total_log[:]
    s.set_prop("log", thislog)
    return s
