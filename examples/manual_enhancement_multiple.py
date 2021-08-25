# {{{ Docstring
"""
Manual Enhancement Multiple
===========================
This script enables the data from multiple spectra to be 
plotted together for comparison under different processing 
conditions (i.e. concentration and heating corrections).

Note: It is recommended that this script only be used
after determining the parameters for each enhancement 
spectrum using the 'manual_enhancement.py' example script,
but can be used to process multiple enhacnements at once ...............

Note: All spectra need to be the same length (i.e. 
len(s.getaxis("power")) in order for this to work). Updates
to get this to work otherwise can be made in the future.

Parameters (in the order they appear)
-------------------------------------
file_location: str - exp_type of your enhancement data file from 
    the _pyspecdata or .pyspecdata file (typically 'ODNP_NMR_comp/ODNP'
    if connected to the G-drive; in this example it is not).

plotname: str - name of the final plot(s) on which all 
    enhancements or k_sigma*s(p), etc. are plotted.

names: list of strings - names for the data which will appear
    in the figure legend. Be sure to put these in order of the
    files you load in later on.

colors: list of strings - colors for the data from the matplotlib xkcd
    library. It's optional to change this.

show_proc: Boolean - determines whether the figures from processing 
    each of the enhancements appears ('True') or just the figures from 
    the last enhancement ('False'). 

plot_all: Boolean - nddata sets will be filled for enhancements (and 
    k_sigma*s(p) curves if applicable) and all of the datasets will be 
    plotted on the same figure with plotname as the title and names[i]
    as the respective legend names.

save_figs: Boolean - determines whether the figures from plot_all will be 
    saved.

savedir: str - directory in which you would like your figures to be saved 
    (make sure there is no "/" on the end).You can put "." for the current 
    working directory if desired, but it's recommended that you have a 
    separate output file directory.

use_T1s: Boolean - determines whether T1 values and corresponding powers
    will be fit to a first order linear function and divided out of the 
    enhancement as a heating correction. If this is True, the function of 
    this script that divides out T1(p=0) if given will be overridden by
    the T1(p) correction. 

auto_T1s: Boolean - determines whether the nddatas that fill a list, are...
    True -- automatically constructed via loading from a csv (columns: 
    "power (W)", "$T_1$ (s)") using pandas.
    Note: There is a script in progress for saving multiple T1 values 
    into a csv automatically, which is where $T_1$ comes from.
    OR
    False --  manually constructed with entered T1 values 
    and powers are used (False is applied for this case, see T1s_x).

T1_filelist: list of strings - list of csv filenames containing the T1 values
    and corresponding power values, no extension. This file should be located 
    in the same 'savedir' directory where output files are allocated to, since
    at some point, this will be a type of output file from IR processing.

T1s_x: nddata - ndshape of [('power',n)], contains the T1 values 
    in s as data and the corresponding n power values in W at which the IR
    experiments were performed as a label.You need a separate T1s_x for each
    enhancement dataset. 

idx: int - index of the dataset starting at 0 and counting up for each
    data set you have. Example: if you have 4 data sets, your indices 
    would be 0, 1, 2, 3. 

filename: str - name of the file containing your enhancement data. No
    extension required.

nodename: str - name of the node in the file containing your data.
    Examples include: "enhancement" "enhancement_1" "signal"

f_range: tuple of floats - frequency range in Hz in which you will 
    crop your data to remove noise. Example: (300,250) #Hz

C: float - concentration in M. If given, you will get (1-E(p))/C. Can
    enter 'None'. Example: 340e-6 #M

T1_0: float - T1(p=0) in s. If given, this script will get (1-E(p))/(C T1(p=0)).
    Can enter 'None' if you don't desire this correction.  You need to enter either
    a value for all enhancement or 'None' for all enhancements. DONT MIX AND MATCH
    because this will result in data that is not directly comparable.
    Example: 0.907 #s

ppt: float - the ratio of f_NMR to f_ESR for your specific type of sample (spin 
    label) in MHz/GHz. Can enter 'None'. Example: 1.5154 #MHz/GHz for MTSL in pR

"""
# }}}
# {{{ Imports & Initialization
from pyspecdata import *
from pyspecProcScripts import *
from sympy import symbols, Symbol, latex, limit, init_printing
from numpy import *
import pandas as pd
import time
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.0),
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
    }
)
logger = init_logging("info")
t2 = symbols("t2")
fl = fl_mod()
measured_vs_actual = 22.0  # how many dB down the split + measured power is from
# the forward power
# }}}
# Input the Parameters
date = time.strftime("%y%m%d")
file_location = "odnp"
postproc = "spincore_ODNP_v2"
plotname = "Q183R1a pR All Environment Comparison"
names = ["DDM KCl", "KI", "$KH_{2}PO_{4}$", "DHPC KCl"]
colors = [
    "dark magenta",
    "bright red",
    "gold",
    "leaf green",
    "cerulean",
    "grape",
    "wine",
]
show_proc = False
plot_all = True
save_figs = True
savedir = "."
if save_figs:
    os.chdir(savedir)
# Load T1_values
use_T1s = True
auto_T1s = False
if auto_T1s:
    T1_filelist = ["210813_Q183R1a_pR_DDM_T1s","210824_Q183R1a_pR_KI_T1s",
            "210824_Q183R1a_pR_KH2PO4_T1s","210824_Q183R1a_pR_DHPC_T1s"]
    T1s = []
    for T1_file in T1_filelist:
        T1_df = pd.read_csv("%s/%s.csv"%(savedir,T1_file))
        T1_data = nddata(T1_df["$T_1$ (s)"].to_numpy(),"power").labels(
                "power", T1_df["power (W)"].to_numpy())
        T1s.append(T1_data)
else:
    T1s_a = nddata(r_[0.959, 1.150, 1.305, 1.394, 1.542], "power").labels(
            "power", r_[0.0, 2.0, 2.51, 3.16, 3.98])
    T1s_b = nddata(r_[1.076, 1.302, 1.401, 1.454, 1.543], "power").labels(
            "power", r_[0.0, 2.0, 2.51, 3.16, 3.98])
    T1s_c = nddata(r_[1.126, 1.461, 1.499, 1.592], "power").labels(
            "power", r_[0.0, 2.51, 3.16, 3.98])
    T1s_d = nddata(r_[1.132, 1.432, 1.464, 1.534, 1.586], "power").labels(
            "power", r_[0.0, 2.0, 2.51, 3.16, 3.98])
    T1s = [T1s_a, T1s_b, T1s_c, T1s_d]

for (idx, filename, nodename, f_range, C, T1_0, ppt) in [
    (
        0,
        "210707_Q183R1a_pR_DDM_ODNP",
        "enhancement",
        (-225, 75),
        207.4e-6,
        0.959,
        1.5154,
    ),
    (
        1,
        "210708_Q183R1a_pR_KI_ODNP",
        "enhancement",
        (-200, 100),
        113.1e-6,
        1.076,
        1.5154,
    ),
    (
        2,
        "210708_Q183R1a_pR_KH2PO4_ODNP",
        "enhancement",
        (-225, 75),
        115.2e-6,
        1.126,
        1.5154,
    ),
    (
        3,
        "210707_Q183R1a_pR_DHPC_ODNP",
        "enhancement",
        (-225, 75),
        103.2e-6,
        1.132,
        1.5154,
    ),
]:
    if use_T1s:
        T1_vals = T1s[idx]
    else: 
        T1_vals = None

    # {{{ Processing Data
    titl = "%s %s" % (" ".join(filename.split("_")[1:-1]), filename.split("_")[0])
    s = find_file(
        filename,
        exp_type=file_location,
        expno=nodename,
        postproc=postproc,
        lookup=lookup_table,
    )
    fl.next("raw data")
    s = s["t2":(-1e3, 1e3)]
    fl.image(s)
    ph0 = s["power", -4].sum("t2")
    ph0 /= abs(ph0)
    s /= ph0
    fl.next("phased")
    fl.image(s)
    s = s["t2":f_range]
    fl.next("sliced to integration bounds")
    fl.image(s)
    if "ph2" in s.dimlabels:
        s = s["ph1", 1]["ph2", 0]
    else:
        s = s["ph1", 1]
    fl.next("sliced to integration bounds\nselect coh. pathway\nreal")
    fl.image(s.real)
    s.integrate("t2")
    s /= max(s.data.real)
    # }}}
    # {{{ Getting Power Axis
    s.setaxis("power", r_[-9999, array(s.get_prop("meter_powers"))])
    print("here are the dBm", s.getaxis("power"))
    s.setaxis("power", lambda x: 1e-3 * 10 ** ((x + measured_vs_actual) / 10.0))
    print("here are the powers", s.getaxis("power"))
    s.set_units("power", "W")

    if plot_all and idx == 0:
        enhancements = nddata(
            zeros((len(names), len(s.getaxis("power"))), dtype="complex128"),
            ["sample", "power"],
        ).labels("power", s.getaxis("power"))
        curves = nddata(
            zeros((len(names), len(s.getaxis("power"))), dtype="complex128"),
            ["sample", "power"],
        ).labels("power", s.getaxis("power"))
        ksigmasp = nddata(
            zeros((len(names), len(s.getaxis("power"))), dtype="complex128"),
            ["sample", "power"],
        ).labels("power", s.getaxis("power"))

    # }}}
    # {{{ Integral
    fl.next("simple integral")
    fl.plot(
        s["power", :-3], "ko", human_units=False
    )  # human_units = False because the range for the red points tries to force mW, which is incompatible with W
    fl.plot(s["power", -3:], "ro", human_units=False)
    plt.ylabel("integrated $^1$H NMR signal")
    # }}}
    # {{{ Enhancement
    s /= s["power", 0]
    fl.next("enhancement curve")
    fl.plot(s["power", :-3], "ko", human_units=False)
    fl.plot(s["power", -3:], "ro", human_units=False)
    plt.xlabel("power (W)")
    plt.ylabel("enhancement")
    if plot_all:
        try:
            enhancements["sample", idx] = s
        except:
            print(
                "The enhancement data you are trying to add to the nddata (idx=%0.0f) has %0.0f power points while the nddata holds exactly%0.0f power points."
                % (idx, len(s.getaxis("power")), len(enhancements.getaxis("power")))
            )

    if use_T1s:
        T1_vals = T1s[idx]
    else: 
        T1_vals = None

    if C is not None:
        epsilon = 1 - s.C
        ks = epsilon / C
        if T1_0 is None:
            fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$" % titl)
        else:
            if T1_vals is None:  # if T1_0 is not None and T1_vals is None
                ks /= T1_0
                fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$" % titl)
            else:
                temp = ks.C / T1_0
                if ppt is not None:
                    temp *= ppt * 1e-3 
                if plot_all: 
                    try:
                        curves["sample", idx] = temp
                    except:
                        print(
                    "The enhancement data you are trying to add to the nddata (idx=%0.0f) has %0.0f power points while the nddata holds exactly%0.0f power  points."
                    % (idx, len(s.getaxis("power")), len(enhancements.getaxis("power")))
                )
                m, b = np.polyfit(T1_vals.getaxis("power"), T1_vals.data, 1)
                T1_p = lambda p: (m * p) + b
                T1_p = nddata(T1_p(ks.getaxis("power")),"power").labels(
                    "power", ks.getaxis("power")
                )
                ks /= T1_p
                fl.next("%s\n$k_{\sigma}s(p)$" % titl)
        if ppt is not None:
            ks *= ppt * 1e-3  # gets you in units of k_sigma!!!
        fl.plot(ks["power", :-3], "ko", human_units=False)
        fl.plot(ks["power", -3:], "ro", human_units=False)
        plt.xlabel("power (W)")
        if T1_0 is None:
            plt.ylabel("$k_{\sigma}s(p)T_{1}(p)$")
        else:
            if T1_vals is None:
                plt.ylabel("$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$")
            else:
                plt.ylabel("$k_{\sigma}s(p)$")
        if plot_all:
            try:
                ksigmasp["sample", idx] = ks
            except:
                print(
                    "The enhancement data you are trying to add to the nddata (idx=%0.0f) has %0.0f power points while the nddata holds exactly%0.0f power  points."
                    % (idx, len(s.getaxis("power")), len(enhancements.getaxis("power")))
                )
    if show_proc:
        fl.show()
# }}}
fl.show()
# {{{ Plotting all results together
if plot_all:
    figure(figsize=(7, 5))
    title("%s %s\nenhancements" % (filename.split("_")[0], plotname))
    enhancements.labels("power", s.getaxis("power"))
    for i, name in enumerate(names):
        plot(
            enhancements["sample", i]["power", :-3],
            color="xkcd:%s" % colors[i],
            marker="o",
            ls="",
            human_units=False,
            label="%s" % name,
            alpha=0.75,
        )
        plot(
            enhancements["sample", i]["power", -3:],
            color="xkcd:%s" % colors[i],
            marker="x",
            ls="",
            human_units=False,
            markersize=7,
        )  # label='%s back'%name,
    plt.ylabel("E")
    plt.xlabel("power (W)")
    plt.legend()  # bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig(
            "%s_%s_enhancement.png" % (date, plotname),
            transparent=True,
            overwrite=True,
            bbox_inches="tight",
        )
    show()
if plot_all and C is not None:
    figure(figsize=(7, 5))
    title("%s %s\n$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(p=0)$" % (
        filename.split("_")[0], plotname))
    curves.labels("power", s.getaxis("power"))
    for i, name in enumerate(names):
        plot(
            curves["sample", i]["power", :-3],
            color="xkcd:%s" % colors[i],
            marker="o",
            ls="",
            human_units=False,
            label="%s" % name,
            alpha=0.75,
        )
        plot(
            curves["sample", i]["power", -3:],
            color="xkcd:%s" % colors[i],
            marker="x",
            ls="",
            human_units=False,
            markersize=7,
        )  # label='%s back'%name,
    plt.ylabel("$k_{\sigma}s(p)$ $\div$ $T_{1}(p=0)$")
    plt.xlabel("power (W)")
    plt.legend()  # bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig(
            "%s_%s_ksgima_sp_T10.png" % (date, plotname),
            transparent=True,
            overwrite=True,
            bbox_inches="tight",
        )
    show()
    figure(figsize=(7, 5))
    title("%s %s\n$k_{\sigma}s(p)$" % (filename.split("_")[0], plotname))
    ksigmasp.labels("power", s.getaxis("power"))
    for i, name in enumerate(names):
        plot(
            ksigmasp["sample", i]["power", :-3],
            color="xkcd:%s" % colors[i],
            marker="o",
            ls="",
            human_units=False,
            label="%s" % name,
            alpha=0.75,
        )
        plot(
            ksigmasp["sample", i]["power", -3:],
            color="xkcd:%s" % colors[i],
            marker="x",
            ls="",
            human_units=False,
            markersize=7,
        )  # label='%s back'%name,
    plt.ylabel("$k_{\sigma}s(p)$ $(s^{-1})$")
    plt.xlabel("power (W)")
    plt.legend()  # bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig(
            "%s_%s_ksgima_sp.png" % (date, plotname),
            transparent=True,
            overwrite=True,
            bbox_inches="tight",
        )

show()
# }}}
