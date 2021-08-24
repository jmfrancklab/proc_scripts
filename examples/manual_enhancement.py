# Docstring
"""
Manual Enhancement
==================
This script allows you to process your enhancement data for 
an individual ODNP experiment manually, in a relatively simple
and shortwinded manner. There are many customizable options
depending on the data you collected and the figures and other
outputs you'd like to save. 

Parameters (in the order they appear)
-------------------------------------
file_location: str - exp_type of your enhancement data file from 
    the _pyspecdata or .pyspecdata file (typically 'ODNP_NMR_comp/ODNP'
    if connected to the G-drive; in this example it is not).

postproc: str - the type of postprocessing you require as written
    in pyspecProcScripts/load_file.py. (This script uses v2 because 
    it accounts for any old 8-step phase-cycled data).

export_csv: Boolean - determines whether a csv is generated from
    the enhancement data (columns: power (W), Re[E(p)], Im[E(p)]).

save_figs: Boolean - determines whether to save enhancement and
    ksigma*s(p) plots. Filename format:
    "{today's date in YYMMDD}_{root filename}_{plot type}.png"
    ex: 
    "210720_150uM_TEMPOL_SMB_enhancement.png"

savedir: str - directory in which you would like your figures to be saved 
    (make sure there is no "/" on the end).You can put "." for the current 
    working directory if desired, but it's recommended that you have a 
    separate output file directory.

auto_T1s: Boolean - determines whether manually entered T1 values 
    and powers are used (False is applied for this case)for the 
    'T1s' nddata: ndshape of [('power', n)] where n is the number
    of T1 values measured OR  if the array is constructed by loading
    from a csv (columns: "power (W)", "$T_1$ (s)") using pandas.
    Note: There is a script in progress for saving multiple T1 values 
    into a csv automatically, which is where $T_1$ comes from.

T1_file: str - name of the csv file containing T1_values, no extension. This
    file should be located in the same 'savedir' directory where output
    files are allocated to, since at some point this will be a type of 
    output file.

manual_T1s: arr - array of T1 values corresponding to power_vals. By
    giving these values, and entering 'T1s' for T1_vals, this script 
    will correct the data by dividing out a first-order fit of T1(p).
    It will do the same if you do auto 

power_vals: arr - array of power values corresponding to T1 values in
    manual_T1s

filename: str - name of the file containing your enhancement data. No
    extension required.

nodename: str - name of the node in the file containing your data.
    Examples include: "enhancement" "enhancement_1" "signal"

f_range: tuple of floats - frequency range in Hz in which you will 
    crop your data to remove noise. Example: (300,250) #Hz

C: float - concentration in M. If given, you will get (1-E(p))/C. Can
    enter 'None'. Example: 340e-6 #M

T1_0: float - T1(p=0) in s. If given, this script will get (1-E(p))/(C T1(p=0)).
    Can enter 'None'. Example: 0.907 #s

T1_vals: nddata - ndshape of [("power", n)] where n is the number of 
    measured T1 values. This is entered as 'T1s' or 'None' depending
    on whethere you have the values and want the correction.

ppt: float - the ratio of f_NMR to f_ESR for your specific type of sample (spin 
    label) in MHz/GHz. Can enter 'None'. Example: 1.5154 #MHz/GHz for MTSL in pR
"""

# {{{ Imports & Initialization
from pyspecdata import *
from pyspecProcScripts import *
from sympy import symbols, Symbol, latex, limit, init_printing
from numpy import *
import pandas as pd
import os
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
# Input parameters & Load Data
date = time.strftime("%y%m%d")
file_location = "odnp"
postproc = "spincore_ODNP_v2"
export_csv = True
save_figs = True
savedir = "C:/Users/saman/Research/Data/pR_ODNP/output_files/"
if save_figs:
    os.chdir(savedir)
auto_T1s = True
if auto_T1s:
    T1_file = "210813_Q183R1a_pR_DDM_T1s"
    T1_df = pd.read_csv("%s/%s.csv" % (savedir, T1_file))
    T1s = nddata(T1_df["$T_1$ (s)"].to_numpy(), "power").labels(
        "power", T1_df["power (W)"].to_numpy()
    )
else:
    manual_T1s = r_[0.959, 1.150, 1.305, 1.394, 1.542]
    power_vals = r_[0.0, 2.00, 2.51, 3.16, 3.98]
    T1s = nddata(manual_T1s, "power").labels("power", power_vals)
# }
for (filename, nodename, f_range, C, T1_0, T1_vals, ppt) in [
    # Example for old 8-step phase cycled data, preprocessed
    # with the spincore_ODNP_v2 processing.
    #    (
    #        "210507_TEMPOL_150uM__cap_probe_DNP_1",
    #        "signal",
    #        (-150, 250),
    #        150e-6,
    #        1.87,
    #        None,
    #        1.5163,
    #    )
    # Example with realistic data with T1 values; you can enter
    # manually above, or import from a csv (ask Sam). You can also
    # choose to enter C, T1_0, T1_vals and/or ppt as None.
    (
        "210707_Q183R1a_pR_DDM_ODNP",
        "enhancement",
        (-225, 75),
        207.4e-6,
        0.959,
        T1s,
        1.5154,
    ),
]:
    # {{{ Load and plot raw data
    outname = "%s_%s" % (date, "_".join(filename.split("_")[1:-1]))
    titl = "%s %s" % (" ".join(filename.split("_")[1:-1]), filename.split("_")[0])
    s = find_file(
        filename,
        exp_type=file_location,
        expno=nodename,
        postproc=postproc,
        lookup=lookup_table,
    )
    fl.next("%s\nraw data" % titl)
    # }}}
    # {{{ Cut off trailing data and phase
    s = s["t2":(-1e3, 1e3)]
    fl.image(s)
    ph0 = s["power", -4].sum("t2")
    ph0 /= abs(ph0)
    s /= ph0
    fl.next("%s\nphased" % titl)
    fl.image(s)
    # }}}
    # {{{ Slice to integration bounds
    s = s["t2":f_range]
    fl.next("%s\nsliced to integration bounds" % titl)
    fl.image(s)
    # }}}
    # {{{ Select coherence pathway
    if "ph2" in s.dimlabels:
        s = s["ph1", 1]["ph2", 0]
    else:
        s = s["ph1", 1]
    fl.next("%s\nselect coh. pathway - real" % titl)
    fl.image(s.real)
    # }}}
    # {{{ Integrate and normalize --> Enhancement
    s.integrate("t2")
    s /= max(s.data.real)
    # }}}
    # {{{ Set power axis
    s.setaxis("power", r_[-9999, array(s.get_prop("meter_powers"))])
    print("here are the dBm", s.getaxis("power"))
    s.setaxis("power", lambda x: 1e-3 * 10 ** ((x + measured_vs_actual) / 10.0))
    print("here are the powers", s.getaxis("power"))
    s.set_units("power", "W")
    # }}}
    # {{{ Simple integral
    fl.next("%s\nsimple integral" % titl)
    fl.plot(s["power", :-3], "ko", human_units=False)
    # (human_units = False because the range for the red points tries to force mW, which is incompatible with W)
    fl.plot(s["power", -3:], "ro", human_units=False)
    plt.xlabel("power (W)")
    plt.ylabel("integrated $^1$H NMR signal")
    # }}}
    # {{{ Enhancement
    s /= s["power", 0]
    fl.next("%s\nenhancement curve" % titl)
    fl.plot(s["power", :-3], "ko", human_units=False)
    fl.plot(s["power", -3:], "ro", human_units=False)
    plt.xlabel("power (W)")
    plt.ylabel("enhancement")
    if save_figs:
        plt.savefig(
            "%s_enhancement.png" % outname,
            transparent=True,
            overwrite=True,
            bbox_inches="tight",
        )
    # }}}
    # {{{ Export enhancement data as a csv (if prompted)
    if export_csv:
        import pandas as pd

        reenhancementdf = pd.DataFrame(data=s.real.data, columns=["Re[E(p)]"])
        imenhancementdf = pd.DataFrame(data=s.imag.data, columns=["Im[E(p)]"])
        enhancementdf = pd.DataFrame(data=s.data, columns=["enhancement"])
        powerdf = pd.DataFrame(data=s.getaxis("power"), columns=["power"])
        enhancement = powerdf.join(reenhancementdf)
        enhancement2 = enhancement.join(imenhancementdf)
        enhancement3 = enhancement2.join(enhancementdf)
        enhancement3.to_csv("%s_enhancement.csv" % outname, index=False)
    # }}}
    # {{{ Correcting the enhancement to obtain ksigma*s(p), etc...
    # Concentration...
    if C is not None:
        epsilon = 1 - s.C
        ks = epsilon / C
        # T1(0)...
        if T1_0 is None:
            fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$" % titl)
        else:
            if T1_vals is None:  # if T1_0 is not None and T1_vals is None
                ks /= T1_0
                fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$" % titl)
            #  T1(p)...
            else:  # (if T1_vals is not None -- we dont want T1_vals to be divided out along with T1_0)
                temp = ks.C / T1_0
                fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$" % titl)
                fl.plot(temp["power", :-3], "ko", human_units=False)
                fl.plot(temp["power", -3:], "ro", human_units=False)
                plt.xlabel("power (W)")
                plt.ylabel("$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$")
                m, b = np.polyfit(T1_vals.getaxis("power"), T1_vals.data, 1)
                T1_p = lambda p: (m * p) + b
                T1_p = nddata(T1_p(ks.getaxis("power")), "power").labels(
                    "power", ks.getaxis("power")
                )
                ks /= T1_p
                fl.next("%s\n$k_{\sigma}s(p)$" % titl)
        # ppt correction for increased sigfigs...
        if ppt is not None:
            ks *= ppt * 1e-3  # (gets you in units of k_sigma)
        # plotting corrected enhancement / ksigma*s(p)...
        fl.plot(ks["power", :-3], "ko", human_units=False)
        fl.plot(ks["power", -3:], "ro", human_units=False)
        plt.xlabel("power (W)")
        # determining the y-label based on which corrections were made...
        if T1_0 is None:
            plt.ylabel("$k_{\sigma}s(p)T_{1}(p)$")
        else:
            if T1_vals is None:
                plt.ylabel("$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$")
            else:
                plt.ylabel("$k_{\sigma}s(p)$")
        # saving the figure if prompted to...
        if save_figs:
            plt.savefig(
                "%s_ksigma.png" % outname,
                transparent=True,
                overwrite=True,
                bbox_inches="tight",
            )
    # }}}
fl.show()
