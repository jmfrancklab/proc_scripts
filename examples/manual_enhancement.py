# {{{ Docstring
'''
Manual Enhancement
==================
This script allows you to process your enhancement data for 
an individual ODNP experiment manually, in a relatively simple
and shortwinded manner. There are many customizable options
depending on the data you collected and the figures and other
outputs you'd like to save. 

Parameters (in the order they appear)
----------
file_location: str - exp_type of your enhancement data file from 
    the _pyspecdata or .pyspecdata file
postproc: str - the type of postprocessing you require as written
    in pyspecProcScripts/load_file.py
export_csv: Boolean - determines whether a csv is generated from
    the enhancement data
auto_T1s: Boolean - determines whether the T1 values and corresponding
    power values are loaded from a csv into an nddata [('power',n)]
    where n is the number of T1 values measured 
    OR if you have to manually enter your T1 and power values.
save_figs: Boolean - determines whether to save enhancement and
    ksigma*s(p) plots.
    Filename format = today's date, root filename, plot type.
        ex: "210720_150uM_TEMPOL_SMB_enhancement.png"
save_location: str - path to where you want your figures saved.
manual_T1s: arr - array of T1 values corresponding to power_vals. This
    script will correct the data by dividing out T1(p) if you give OR
    load T1_vals.
power_vals: arr - array of power values corresponding to T1 values in
    manual_T1s
filename: str - name of the file containing your data
nodename: str - name of the node in the file containing your data
f_range: tuple of floats - frequency range in MHz
C: float - concentration in M. If given, you will get (1-E(p))/C. Can
    enter None.
T1_0: float - T1(p=0) in s. If given, this script will get (1-E(p))/(C T1(p=0)).
    Can enter None.
ppt: float - the ratio of f_NMR to f_ESR for your specific type of sample (spin 
    label) in MHz/GHz. 
'''
# }}}
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
# {{{ Input parameters & Load Data
file_location = "odnp"
postproc = "spincore_ODNP_v2" 
export_csv = True
auto_T1s = True
save_figs = False 
save_location = "C:/Users/saman/Research/Data/pR_ODNP/output_files/"
if save_figs:
    os.chdir(save_location)
date = time.strftime("%y%m%d")
# { Load T1_values
if auto_T1s:
    path = save_location
    T1_df = pd.read_csv('%s210707_Q183R1a_pR_DDM_T1s.csv'%path)
    print(T1_df);quit()
    T1s = T1_df['T1']
else:
    manual_T1s = r_[]
    power_vals = r_[]
    T1s = nddata(manual_T1s,'power').labels('power',power_vals)
# }
for (filename, nodename, f_range, C, T1_0, T1_vals, ppt) in [
    # Example for old 8-step phase cycled data, preprocessed
    # with the spincore_ODNP_v2 processing. 
    (
        "210507_TEMPOL_150uM__cap_probe_DNP_1",
        "signal",
        (-150, 250),
        150e-6,
        1.87,
        None,
        1.5163,
    )
    # Example for data with T1 values, you can enter manually
    # above, or import from a csv. 
    (
        "210714_150uM_TEMPOL_SMB_ODNP",
        "enhancement_real",
        (-200,200),
        146.7e-6,
        3.56,
        T1s,
        1.51563,
    )
    #        ('210707_Q183R1a_pR_DDM_ODNP','enhancement',
    #            (-225,75),207.4e-6,None,None,1.5154), # have T1_values for this data
]:
    # }}}
    outname = "%s_%s" % (date, "_".join(filename.split("_")[1:-1]))
    titl = "%s %s" % (" ".join(filename.split("_")[1:-1]), filename.split("_")[0])
    # {{{ Load and plot raw data
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
    fl.plot(
        s["power", :-3], "ko", human_units=False
    )
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
                m, b = np.polyfit(T1_vals.getaxis("power"), T1_vals.data, 1)
                T1_p = lambda p: (m * p) + b
                T1_p = nddata(T1_p(ks.getaxis("power")), ks.getaxis("power")).labels(
                    "power", ks.getaxis("power"))
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
