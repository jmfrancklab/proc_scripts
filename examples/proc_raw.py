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
import numpy as np
import sys
import os
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)

assert len(sys.argv) == 4
d = find_file(
    sys.argv[2], exp_type=sys.argv[3], expno=sys.argv[1], lookup=lookup_table
)
print("postproc_type:", d.get_prop("postproc_type"))
with figlist_var() as fl:
    d.squeeze()
    print("=" * 13 + "ACQ PARAMS" + "=" * 13)
    for k, v in d.get_prop("acq_params").items():
        print(f"{k:>25s} : {v}")
    fl.next("raw data")
    print("=" * 36)

    def image_or_plot(d):
        if len(d.dimlabels) == 1:
            fl.plot(d)
        elif len(d.dimlabels) == 2:
            iterdim = d.shape.min()
            if d.shape[iterdim] > 5:
                d.pcolor()
                return
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
            rows = np.prod([d.shape[j] for j in d.dimlabels[:-1]])
            if rows < 500:
                fl.image(d)
            else:
                fl.image(d, interpolation="bilinear")

    image_or_plot(d)
    if "nScans" in d.dimlabels:
        d.mean("nScans")
        fl.next("signal averaged along nScans")
        image_or_plot(d)
    if d.get_prop("coherence_pathway") is not None:
        fl.next("sum of abs of all coherence pathways (for comparison)")
        forplot = abs(d)
        guess_direct = (
            d.shape.max()
        )  # guess that the longest dimension is the direct
        if guess_direct == "indirect":
            temp = d.shape
            temp.pop("indirect")
            guess_direct = temp.max()
        forplot.mean_all_but(
            list(d.get_prop("coherence_pathway").keys()) + [guess_direct]
        )
        image_or_plot(forplot)
        d = select_pathway(d, d.get_prop("coherence_pathway"))
        fl.next("with coherence pathway selected")
        image_or_plot(d)
