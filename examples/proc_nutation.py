"""
Process nutation data
====================
`py proc_nutation.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.

Tested with:

``py proc_nutation.py nutation_2 240710_27mM_TEMPOL_chokes_T_atprobe_nutation.h5 ODNP_NMR_comp/nutation``

"""
import pyspecProcScripts as prscr
import pyspecdata as psd
import sys
import os
from itertools import cycle

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)

assert len(sys.argv) == 4
d = find_file(
    sys.argv[2],
    exp_type=sys.argv[3],
    expno=sys.argv[1],
    lookup=prscr.load_data.lookup_table,
)
if d.get_prop("postproc_type") == "spincore_SE_v1":
    d.rename("indirect", "p_90")
    d["p_90"] *= 1e-9
elif "indirect" in d.dimlabels:
    d.rename("indirect", "p_90")
with psd.figlist_var() as fl:
    # {{{ Manual processing
    d.squeeze()
    d.ift("t2")
    d["t2"] -= 2e-3  # eyeball correction
    d.ft("t2")
    fl.next("Raw Data")
    fl.image(d, interpolation="auto")
    fl.next("raw data, phase")
    ph0 = (
        prscr.select_pathway(d)["t2":(-100, 0)]["p_90", :5]
        .sum("t2")
        .sum("p_90")
    )
    ph0 /= abs(ph0)
    fl.image(
        prscr.select_pathway(d / ph0).angle * 180 / pi, interpolation="auto"
    )
    fl.next("raw data, abs")
    prscr.select_pathway(d).C.mean("nScans").run(abs).pcolor()
    fl.next("time domain")
    d.ift("t2")
    fl.image(d["t2":(None, 20e-3)], interpolation="auto")
    d_manualslice = d.C["t2":(0, None)]
    d.ft("t2")
    d_manualslice["t2", 0] *= 0.5
    d_manualslice.ft("t2")
    d_manualslice /= prscr.zeroth_order_ph(
        prscr.select_pathway(d_manualslice)["t2":(-250, 250)].mean("nScans")
    )
    my_signs = prscr.determine_sign(
        prscr.select_pathway(d_manualslice)["t2":(-250, 250)].mean("nScans")
    )
    fl.next("manual fid slice")
    fl.image(
        prscr.select_pathway(d_manualslice)["t2":(-250, 250)],
        interpolation="auto",
    )
    fl.next("manual fid slice, with sign flip")
    fl.image(
        prscr.select_pathway(d_manualslice)["t2":(-250, 250)] * my_signs,
        interpolation="auto",
    )
    # }}}
    ## {{{ correct for evolution during pulse
    ##     I do this, but if I plot it, I do find that it's not significant
    # d *= np.exp(
    #    -1j * 2 * pi * (2 * d.fromaxis("p_90") / pi) * d.fromaxis("t2")
    # )
    ## }}}
    # {{{ do it the "correct" way
    fl.next("prscr.fid_from_echo result")
    d = prscr.fid_from_echo(d)
    d /= prscr.zeroth_order_ph(
        prscr.select_pathway(d)["t2":(-250, 250)].mean("nScans")
    )
    fl.image(prscr.select_pathway(d)["t2":(-250, 250)], interpolation="auto")
    fl.next("prscr.fid_from_echo result, with sign flip")
    my_signs = prscr.determine_sign(
        prscr.select_pathway(d)["t2":(-250, 250)].mean("nScans")
    )
    fl.image(
        prscr.select_pathway(d)["t2":(-250, 250)] * my_signs,
        interpolation="auto",
    )
    # }}}
    if d.get_prop("coherence_pathway") is not None:
        d.mean("nScans")
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
        fl.plot(forplot)
        d = prscr.select_pathway(d)
        fl.next("with coherence pathway selected")
        d["t2":(-300, 300)].pcolor()
