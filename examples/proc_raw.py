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
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)

assert len(sys.argv) == 4
d = find_file(
    sys.argv[2], exp_type=sys.argv[3], expno=sys.argv[1], lookup=lookup_table
)
print("postproc_type:",d.get_prop('postproc_type'))
with figlist_var() as fl:
    d.squeeze()
    fl.next("raw data")

    def image_or_plot(d):
        if len(d.dimlabels) == 1:
            fl.plot(d)
        elif len(d.dimlabels) == 2:
            iterdim = d.shape.min()
            untfy_axis = d.unitify_axis(iterdim)
            for idx in range(d.shape[iterdim]):
                c = next(colorcyc)
                fl.plot(
                    d[iterdim, idx],
                    label=f"{untfy_axis}={d[iterdim][idx]}",
                    c=c,
                    alpha=0.5,
                    human_units=False,
                )
                fl.plot(
                    d[iterdim, idx].imag,
                    label=f"{untfy_axis}={d[iterdim][idx]}",
                    c=c,
                    alpha=0.1,
                    human_units=False,
                )
        else:
            fl.image(d)

    image_or_plot(d)
    if "nScans" in d.dimlabels:
        d.mean("nScans")
        fl.next("signal averaged along nScans")
        image_or_plot(d)
    if d.get_prop("coherence_pathway") is not None:
        fl.next("sum of abs of all coherence pathways (for comparison)")
        forplot = abs(d)
        guess_direct = d.shape.max() # guess that the longest dimension is the direct
        forplot.mean_all_but(list(d.get_prop("coherence_pathway").keys())+[guess_direct])
        image_or_plot(forplot)
        d = select_pathway(d, d.get_prop("coherence_pathway"))
        fl.next("with coherence pathway selected")
        image_or_plot(d)
        if len(d.dimlabels) == 2:
            fl.next("pcolor of selected coherence pathway")
            d.pcolor()
