"""
Show data with postproc
=======================
`py proc_raw_data.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are
performed according to the `postproc_type` attribute of the data node.  This
script plots the result, as well as signal that's averaged along the `nScans`
dimension.

Tested with:

``py proc_raw.py echo_6 240620_200uM_TEMPOL_pm_echo.h5 ODNP_NMR_comp/Echoes``

``py proc_raw.py echo_8 240620_200uM_TEMPOL_pm_generic_echo.h5 \
ODNP_NMR_comp/Echoes``

``py proc_raw.py CPMG_9 240620_200uM_TEMPOL_pm_generic_CPMG.h5 \
ODNP_NMR_comp/Echoes``

``py proc_raw.py field_3 240920_27mM_TEMPOL_debug_field \
ODNP_NMR_comp/field_dependent``

``py proc_raw.py ODNP K42.*A1_kRasbatch240814 ODNP_NMR_comp/ODNP``

``py proc_raw.py FIR_34dBm K42.*A1_kRasbatch240814 ODNP_NMR_comp/ODNP``

"""

import pyspecProcScripts as prscr
import sys, os
import matplotlib.pyplot as plt
from itertools import cycle
import pyspecdata as psd

if (
    "SPHINX_GALLERY_RUNNING" in os.environ
    and os.environ["SPHINX_GALLERY_RUNNING"] == "True"
):
    sys.argv = [
        sys.argv[0],
        "echo_6",
        "240620_200uM_TEMPOL_pm_echo.h5",
        "ODNP_NMR_comp/Echoes",
    ]

psd.init_logging(level="debug")

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)

assert len(sys.argv) == 4
d = psd.find_file(
    sys.argv[2],
    exp_type=sys.argv[3],
    expno=sys.argv[1],
    lookup=prscr.lookup_table,
)
print("postproc_type:", d.get_prop("postproc_type"))


def get_first_field_if_structured(d):
    for j in d.dimlabels:
        arr = d[j]
        if arr is not None:
            if arr.dtype.names:  # True if structured array
                print(
                    "for the '",
                    j,
                    "' dimension, you will want to select the field ",
                    arr.dtype.names,
                    " that you want to use, but I'm just picking the first",
                )
                retval = arr[arr.dtype.names[0]]
                if "time" in arr.dtype.names[0].lower():
                    print(
                        "this is called 'time', so I'm assuming it has units"
                        " of seconds"
                    )
                    d.set_units(j, "s")
                if retval[1] / abs(retval[1] - retval[0]) > 10:
                    # there is a large offset to all the numbers
                    retval -= retval[0]
                    print(
                        "I'm also making this axis relative, because it has a"
                        " large offset (it's probably a time axis)"
                    )
                d[j] = retval
    return d


with psd.figlist_var() as fl:
    d.squeeze()
    print("=" * 13 + "ACQ PARAMS" + "=" * 13)
    for k, v in d.get_prop("acq_params").items():
        print(f"{k:>25s} : {v}")
    print("=" * 36)

    def image_or_plot(d):
        if len(d.dimlabels) == 1:
            fl.plot(d)
        elif len(d.dimlabels) == 2:
            iterdim = d.shape.min()
            if d.shape[iterdim] > 3:
                # so that we can do pcolor, if the indirect is a structured
                # array, just pull the first field
                if d[d.dimlabels[0]].dtype.names is not None:
                    the_field = d[d.dimlabels[0]].dtype.names[0]
                    d[d.dimlabels[0]] = d[d.dimlabels[0]][the_field]
                    if "time" in the_field:
                        d[d.dimlabels[0]] -= d[d.dimlabels[0]][0]
                        d[d.dimlabels[0]] /= 60
                        d.set_units(d.dimlabels[0], "min")
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
            d = get_first_field_if_structured(d)
            fl.DCCT(d)

    fl.next("raw data")
    print("about to image or plot", d.shape)
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
        d = prscr.select_pathway(d, d.get_prop("coherence_pathway"))
        fl.next("with coherence pathway selected")
        image_or_plot(d)
        plt.gcf().suptitle("select " + str(d.get_prop("coherence_pathway")))
