# {{{ Imports & Initialization
from pyspecdata import *
from pyspecProcScripts import *
from sympy import symbols, Symbol, latex, limit, init_printing
from numpy import *
import time
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
        "axes.facecolor": (
            1.0,
            1.0,
            1.0,
            0.0,
        ),  # clear ## 1,1,1,0.9 == 90% transparent white
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
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
postproc = "spincore_ODNP_v2"  # 'spincore_ODNP_v1'
export_csv = True
save_figs = False  # saves enhancement and ksigma plots only, youll have to do the
# rest manually. Format todays date, root filename, figure type:
# ex: "210720_150uM_TEMPOL_SMB_enhancement.png'
if save_figs:
    import os

    os.chdir("C:/Users/saman/Research/Sam_Notebooks/ODNP_proc/")
date = time.strftime("%y%m%d")
# { Load T1_values
# import pandas as pd
# path = 'C:/Users/saman/Research/Sam_Notebooks/ODNP_proc'
# T1_df = pd.read_csv('%s/210707_Q183R1a_pR_DDM_T1s.csv'%path)
# print(T1_df);quit()
# }
for (filename, nodename, f_range, C, T1_0, T1_vals, ppt) in [
    (
        "210507_TEMPOL_150uM__cap_probe_DNP_1",
        "signal",
        (-150, 250),
        150e-6,
        1.87,
        None,
        1.5163,
    )
    #        ('210714_150uM_TEMPOL_SMB_ODNP','enhancement_real',
    #            (-200,200),150e-6,3.56,None,1.5163),
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
    )  # human_units = False because the range for the red points tries to force mW, which is incompatible with W
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
    #    plt.show()
    #    plt.close()
    # }}}
    # {{{ Determining ksigma*s(p) or any other possible corrections...
    # { Concentration
    if C is not None:
        epsilon = 1 - s.C
        ks = epsilon / C
        # }
        # { T1(0)
        if T1_0 is None:
            fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$" % titl)
        else:
            if T1_vals is None:  # if T1_0 is not None and T1_vals is None
                ks /= T1_0
                fl.next("%s\n$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$" % titl)
            # }
            # { T1(p)
            else:  # if T1_vals is not None (we dont want T1_vals to be divided out along with T1_0)
                m, b = np.polyfit(T1_vals.getaxis("power"), T1_vals.data, 1)
                T1_p = lambda p: (m * p) + b
                T1_p = nddata(T1_p(ks.getaxis("power")), ks.getaxis("power")).labels(
                    "power", ks.getaxis("power")
                )
                ks /= T1_p
                fl.next("%s\n$k_{\sigma}s(p)$" % titl)
        # }
        # { ppt correction for increased sigfigs
        if ppt is not None:
            ks *= ppt * 1e-3  # gets you in units of k_sigma!!!
        # }
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
        if save_figs:
            plt.savefig(
                "%s_ksigma.png" % outname,
                transparent=True,
                overwrite=True,
                bbox_inches="tight",
            )
    #        plt.show()
    #        plt.close()
    # }}}
    # }
    # { Exporting data from the enhancement curve as csv
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
# }
fl.show()
