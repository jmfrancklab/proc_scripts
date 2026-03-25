"""
Show data with postproc
=======================
`pyspecProcScripts raw EXP_TYPE FILENAME NODENAME`

Fourier transforms (and any needed data corrections for older data) are
performed according to the `postproc_type` attribute of the data node.  This
script plots the result, as well as signal that's averaged along the `nScans`
dimension.

Tested with:

``pyspecProcScripts raw ODNP_NMR_comp/Echoes \
    240620_200uM_TEMPOL_pm_echo.h5 echo_6``

``pyspecProcScripts raw ODNP_NMR_comp/Echoes \
240620_200uM_TEMPOL_pm_generic_echo.h5 echo_8``

``pyspecProcScripts raw ODNP_NMR_comp/Echoes \
240620_200uM_TEMPOL_pm_generic_CPMG.h5 CPMG_9``

``pyspecProcScripts raw ODNP_NMR_comp/field_dependent \
240920_27mM_TEMPOL_debug_field field_3``

``pyspecProcScripts raw ODNP_NMR_comp/ODNP \
K42.*A1_kRasbatch240814 ODNP``

``pyspecProcScripts raw ODNP_NMR_comp/ODNP \
K42.*A1_kRasbatch240814 FIR_34dBm``
"""

from itertools import cycle

import matplotlib.pyplot as plt
import pyspecdata as psd

import pyspecProcScripts as prscr


def image_or_plot(d, fl, colorcyc):
    if len(d.dimlabels) == 1:
        fl.plot(d)
    elif len(d.dimlabels) == 2:
        iterdim = d.shape.min()
        if d.shape[iterdim] > 3:
            # so that we can do pcolor, if the indirect is a structured array,
            # just pull the first field
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
        for dim_name in d.dimlabels:
            arr = d[dim_name]
            if arr is not None and arr.dtype.names:
                print(
                    "for the '",
                    dim_name,
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
                    d.set_units(dim_name, "s")
                if retval[1] / abs(retval[1] - retval[0]) > 10:
                    # there is a large offset to all the numbers
                    retval -= retval[0]
                    print(
                        "I'm also making this axis relative, because it has a"
                        " large offset (it's probably a time axis)"
                    )
                d[dim_name] = retval
        fl.DCCT(d)


def run_raw(exp_type, filename, node):
    psd.init_logging(level="debug")
    colorcyc = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    d = psd.find_file(
        filename,
        exp_type=exp_type,
        expno=node,
        lookup=prscr.lookup_table,
    )
    print("postproc_type:", d.get_prop("postproc_type"))
    with psd.figlist_var() as fl:
        d.squeeze()
        print("=" * 13 + "ACQ PARAMS" + "=" * 13)
        for k, v in d.get_prop("acq_params").items():
            print(f"{k:>25s} : {v}")
        print("=" * 36)
        fl.next("raw data")
        print("about to image or plot", d.shape)
        image_or_plot(d, fl, colorcyc)
        if "nScans" in d.dimlabels:
            d.mean("nScans")
            fl.next("signal averaged along nScans")
            image_or_plot(d, fl, colorcyc)
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
            image_or_plot(forplot, fl, colorcyc)
            d = prscr.select_pathway(d, d.get_prop("coherence_pathway"))
            fl.next("with coherence pathway selected")
            image_or_plot(d, fl, colorcyc)
            plt.gcf().suptitle(
                "select " + str(d.get_prop("coherence_pathway"))
            )


def main(argv=None):
    if argv is None:
        import sys

        argv = sys.argv[1:]
    assert len(argv) == 3
    return run_raw(argv[0], argv[1], argv[2])
