"""
Process nutation data
====================
`py proc_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.

Tested with:

``py proc_nutation.py nutation_2 240710_27mM_TEMPOL_chokes_T_atprobe_nutation.h5 ODNP_NMR_comp/nutation``

"""
from pyspecProcScripts.load_data import lookup_table
from pyspecProcScripts import select_pathway, determine_sign
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
if d.get_prop("postproc_type") == "spincore_SE_v1":
    d.rename("indirect", "p_90")
    d["p_90"] *= 1e-9
with figlist_var() as fl:
    d.squeeze()
    fl.next("raw data")

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
        fl.next("time domain")
        d.ift("t2")
        image_or_plot(d["t2":(None, 20e-3)])
        d.ft("t2")
        # {{{ correct for evolution during pulse
        #     I do this, but if I plot it, I do find that it's not significant
        d *= np.exp(
            -1j * 2 * pi * (2 * d.fromaxis("p_90") / pi) * d.fromaxis("t2")
        )
        # }}}
        my_signs = determine_sign(select_pathway(d)["t2":(-250, 250)])
        d *= my_signs
        fl.next("frequency domain -- after correct and sign flip")
        image_or_plot(d["t2":(-250, 250)])
        fl.next("time domain -- after correct and sign flip")
        d.ift("t2")
        image_or_plot(d["t2":(None, 20e-3)])
        d.ft("t2")
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
        d = select_pathway(d)
        fl.next("with coherence pathway selected")
        image_or_plot(d)
