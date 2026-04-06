import pyspecdata as psd
from pyspecdata import gammabar_H
from .convert_to_power import load_log_data
from .load_data import lookup_table
from .simple_functions import select_pathway


def plot_field(
    data,
    filename,
    exp_type,
    nodename,
    fl=None,
    node_name="log",
):
    log_array = load_log_data(
        filename,
        exp_type,
        nodename=nodename,
        node_name=node_name,
        lookup=lookup_table,
    )
    if "field" not in log_array.dtype.names:
        print("No 'field' column in log data, cannot plot field dependence")
        return
    data = select_pathway(data, data.get_prop("coherence_pathway"))
    carrier = data.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
    data.setaxis(
        "indirect", log_array["field"] * gammabar_H - carrier
    ).set_units("indirect", "Hz")
    forplot = data.ft("t2").run(abs).max("t2")
    if fl is not None:
        fl.next("Peak Offset vs Magnetic Field")
        fl.plot(forplot, "o")
    return forplot
