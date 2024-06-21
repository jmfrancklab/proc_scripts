"""
Show data with postproc
====================
`py proc_raw_data.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.

Tested with:

``py proc_raw.py echo_6 240620_200uM_TEMPOL_pm_echo.h5 ODNP_NMR_comp/Echoes``

``py proc_raw.py echo_8 240620_200uM_TEMPOL_pm_generic_echo.h5 ODNP_NMR_comp/Echoes``

``py proc_raw.py CPMG_9 240620_200uM_TEMPOL_pm_generic_CPMG.h5 ODNP_NMR_comp/Echoes``

"""
from pyspecProcScripts.load_data import lookup_table
from pyspecProcScripts import select_pathway
from pyspecdata import *
import sys
import os

assert len(sys.argv) == 4
d = find_file(sys.argv[2], exp_type=sys.argv[3], expno=sys.argv[1], lookup=lookup_table)
with figlist_var() as fl:
    d.squeeze()
    fl.next("raw data")

    def image_or_plot(d):
        if len(d.dimlabels) > 2:
            fl.image(d)
        else:
            fl.plot(d, alpha=0.5)
            fl.plot(d.imag, alpha=0.1)

    image_or_plot(d)
    if "nScans" in d.dimlabels:
        d.mean("nScans")
        fl.next("signal averaged along nScans")
        image_or_plot(d)
    if d.get_prop("coherence_pathway") is not None:
        forplot = abs(d)
        for thisdim in d.get_prop("coherence_pathway").keys():
            forplot.sum(thisdim)
        fl.next("sum of abs of all coherence pathways (for comparison)")
        image_or_plot(forplot)
        d = select_pathway(d, d.get_prop("coherence_pathway"))
        fl.next("with coherence pathway selected")
        image_or_plot(d)
