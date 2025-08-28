"""
Show data collected from an ODNP experiment (Enhancement and IR) using either
rough or correlated alignment processing functions
==================================================
'py proc_alignedtoi.py FIR_35dBm 250321_280uM_OHT_ODNP_1 ODNP_NMR_comp/ODNP correlation'

Meant to be run from the command line similar to proc_raw and other scripts in
proc_scripts/examples to quickly process and view data collected during an ODNP
experiment

frq_mask() and coherence_unmask_fn() functions taken from correlation_alignment_example.py

WK, 8/13/25 - want to totally gut and rewrite after convo with JF. There are
things to keep at the start with calling a file in the bash terminal. Now want
to edit more of the newly named 'table_of_integrals' function adding an
attribute to choose which type of alignment is being applied but otherwise treat data the same and have a similar resulting plot.  


"""

import pyspecProcScripts as prscr
import pyspecdata as psd
import sys, os, datetime
import matplotlib.pyplot as plt
from itertools import cycle
import numpy as np

psd.init_logging(level="debug")

#{{{ Plotting setup
plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2
plt.rcParams.update({
    "errorbar.capsize": 2,
    "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor": (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
    "figure.figsize": (6, 5),
})

colorcyc_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colorcyc = cycle(colorcyc_list)
# }}}

# {{{ Mask funcs from correl align example
# These (mask and unmask funcs) may need to be moved to within
# 'table_of_integrals' func or as actual funcs for used in and outside of the
# 'toi' function
def frq_mask(s, sigma=150.0):
    """Note that we assume that our mask is a product of a
    frequency-domain and a coherence-domain function.
    This returns a copy multiplied by the square root of the
    frequency-domain part,
    leaving the original data untouched.

    Parameters
    ==========
    s : nddata
        Signal, given in the frequency domain and coherence transfer
        (*vs.* phase) domain.
        The property `coherence_pathway` must be set.
    """
    assert s.get_ft_prop("t2")
    assert s.get_ft_prop(list(s.get_prop("coherence_pathway").keys())[0])
    # {{{ find center frequency
    nu_center = (
        prscr.select_pathway(s, s.get_prop("coherence_pathway"))
        .mean("repeats")
        .argmax("t2")
    )
    # }}}
    # {{{ Make mask using the center frequency and sigma.
    #     Standard gaussian is 2σ² in the denominator -- the extra 2 is
    #     for sqrt.
    frq_mask = np.exp(-((s.fromaxis("t2") - nu_center) ** 2) / (4 * sigma**2))
    # }}}
    # note that when we multiply, we automatically generate a copy
    return s * frq_mask


def coherence_unmask_fn(coh_array):
    """Filters out all but the signal pathway and the {"ph1":0} or
    {"ph1":0,"ph2":0} pathways (depending on which experiment below is used).
    Note this serves as an example function and other filter functions could
    alternatively be used"""

    def set_pathway_true(pathway_dict):
        for j, (k, v) in enumerate(pathway_dict.items()):
            # the last element needs to be treated differently
            if j < len(pathway_dict) - 1:
                thisslice = coh_array[k, v]
            elif len(pathway_dict) == 1:
                coh_array[k, v] = 1
            else:
                thisslice[k, v] = 1

    set_pathway_true(coh_array.get_prop("coherence_pathway"))
    set_pathway_true(
        {k: 0 for k in coh_array.get_prop("coherence_pathway").keys()}
    )
    return coh_array
# }}}

# {{{ Most likely want to keep...
assert len(sys.argv) == 5
filename=sys.argv[2]
exptype=sys.argv[3]
nodename=sys.argv[1]
alignment=sys.argv[4]

dataset = psd.find_file(
    filename,
    exp_type=exptype,
    expno=nodename,
    lookup=prscr.lookup_table,
)

with psd.figlist_var() as fl:
    s = dataset.C
    s = s.squeeze()
    print(s)
    fl.next("raw data")
    fl.image(s)
    signal_range = None
    signal_pathway = None
    direct="t2"
    if nodename == "ODNP": # may be a better way to do this
        indirect="indirect"
    if nodename.startswith("FIR_"):    
        indirect="vd"
    expansion=2
    peak_lower_thresh=0.1
    if signal_range is None:
        center_of_range, half_range = prscr.find_peakrange(
            s, fl=fl, direct=direct, peak_lower_thresh=peak_lower_thresh
        )
        signal_range = s.get_prop("peakrange")
    else:
        if signal_range == "peakrange":
            signal_range = s.get_prop("peakrange")
        center_of_range = np.mean(signal_range)
        half_range = 0.5 * np.diff(signal_range).item()
    signal_range_expanded = (
        center_of_range + expansion * np.r_[-1, 1] * half_range
    )
    signal_pathway = signal_pathway or s.get_prop("coherence_pathway")
    if nodename == "ODNP":
        orig_axis = s["indirect"] # Save the original indirect axis with start and
        #                           stop times for convert_to_power
        s["indirect"] = (
            s["indirect"]["start_times"] - s["indirect"]["start_times"][0]
        )
        s.set_units("indirect", "s")
    print("s before split in processing", s.shape)
    # {{{ Rough Alignment
    if alignment == "rough":
        clock_correction=False # I'm not sure if this is needed
        s.set_plot_color(
            "g"
        )  # this affects the 1D plots, but not the images, etc.
        if nodename.startswith("FIR_"):
            print("I am seeing nodename as 'FIR_*'")
            if clock_correction:
                s = prscr.clock_correct(s)
        fl.next("rough aligned")
        s, _ = prscr.table_of_integrals(s, fl=fl, title="rough")
        fl.next("after RTOI %s" % nodename)
        fl.plot(s)
        if nodename == "ODNP":
            s.set_error(s["indirect", 0].item() * 0.01)  # Not calculating
            #                                              error, just looking
            #                                              for relative errors
            s /= s["indirect", 0:1]
            fl.next("normalized $E(p(t))$ %s" % filename)
            print("s[indirect][0] is ", s["indirect"][0])
            print("s[indirect][1] is ", s["indirect"][1])
            fl.plot(s, "o")
            # {{{ this is just matplotlib time formatting
            ax = plt.gca()
            ax.xaxis.set_major_formatter(
                plt.FuncFormatter(lambda x, _: str(datetime.timedelta(seconds=x)))
            )
            # }}}
            s.set_prop("filename", filename)
            thisfile = s.get_prop("filename")
            print("before convert to P", s)
            s["indirect"] = orig_axis # Reset "indirect" axis back to start stop times
            s = prscr.convert_to_power(s, thisfile, exptype, fl=fl)
            fl.next("normalized $E(p)$ %s" % filename)
            fl.plot(s, "o")
    # }}}
    # {{{ Correlation Aligned
    if alignment == "correlation":
        print("nodename = ", nodename, "\n nodename type is ", type(nodename))
        print("s['t2'] is ", s["t2"])
        s.ift("t2")
        s /= prscr.zeroth_order_ph(
            prscr.select_pathway(s, signal_pathway)
        )
        print("s looks like ", s)
        s["t2"] -= s.getaxis("t2")[0]  # needed for Hermitian Function
        fl.next("first step corr data")
        fl.DCCT(s)
        if nodename == "ODNP":
            best_shift = prscr.hermitian_function_test(
                prscr.select_pathway(s.C.mean("indirect"), signal_pathway)
            ) 
        if nodename.startswith("FIR_"):
            best_shift = prscr.hermitian_function_test(
                prscr.select_pathway(s.C.mean("vd"), signal_pathway)
            )    
        if "nScans" in s.dimlabels:
            print("I am in nScans")
            s.mean("nScans")
            fl.next("signal averaged along nScans")
            fl.plot(s)
        if s.get_prop("coherence_pathway") is not None:
            print("I am in coherence pathways")
            fl.next("sum of abs of all coherence pathways (for comparison)")
            forplot = abs(s)
            guess_direct = (
                s.shape.max()
            )  # guess that the longest dimension is the direct
            if guess_direct == "indirect":
                temp = s.shape
                temp.pop("indirect")
                guess_direct = temp.max()
            forplot.mean_all_but(
                list(s.get_prop("coherence_pathway").keys()) + [guess_direct]
            )
        print(best_shift)
        s.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        s.ft("t2")
        fl.next("Phased and Centered, %s" % nodename)
        fl.DCCT(s)
        mysgn = (
            prscr.select_pathway(s, signal_pathway)
            .C.real.sum("t2")
            .run(np.sign)
        )
        fl.next("Apply correl_align")
        opt_shift = prscr.correl_align(
            s * mysgn,
            frq_mask_fn=frq_mask,
            coherence_unmask_fn=coherence_unmask_fn,
            repeat_dims=indirect,
            max_shift=300, # 3kHz Gaussian mask
            fl=fl,
        )    
        s.ift("t2").ift(list(signal_pathway))
        s *= np.exp(-1j * 2 * np.pi * opt_shift * s.fromaxis("t2"))
        s.ft("t2").ft(
            list(signal_pathway.keys())
        )
        fl.next("Aligned data ($\\nu$)")
        fl.DCCT(s)
        s.ift("t2")
        fl.next("Aligned data ($t$)")
        fl.DCCT(s)

fl.show()
