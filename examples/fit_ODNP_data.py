"""Fitting ODNP Datasets for Ksigma
===================================
The T1(p) and E(p) integrals are pulled from an H5 file and are fit to extract the cross relaxivity of the sample. The R1(p) is first interpolated by fitting the self relaxivity to a polynomial. This fit is then used to fit the linear regime of the enhancement data. The cross relaxivity is then calculated using the output of these fits as well as the sample parameters that are fed to it (e.g., ppt value, and concentration).
"""
from pyspecdata import *
from sympy import symbols, Symbol, latex

target_directory = os.path.normpath(getDATADIR("AG_processed_data"))
h5_folder = "ras.h5"
filename = "220616_E37"
# {{{Dataset Parameters
SL_conc_M = 455.2223e-6
ppt = 1.5167e-3
Ep_pts = 18
# }}}
# {{{ plotting fn
def list_symbs(f):
    # {{{ this is just to show all the parameters
    list_symbs = []
    for j, k in f.output().items():
        s_repr = latex(Symbol(j))
        list_symbs.append(f"${s_repr} = {k:0.5g}$")
    list_symbs = "\n".join(list_symbs)
    # }}}
    return list_symbs


# }}}
with figlist_var() as fl:
    # {{{ load data
    Ep = nddata_hdf5(f"{h5_folder}/{filename}/Ep", directory=target_directory)
    T1p = nddata_hdf5(f"{h5_folder}/{filename}/T1p", directory=target_directory)
    # }}}
    # {{{Plot Ep
    fl.next("E(p)")
    fl.plot(Ep["power", :-3], "o", label="Experimental Data", capsize=6, alpha=0.5)
    fl.plot(Ep["power", -3:], "rx", label="Returning Power Check", capsize=6, alpha=0.5)
    plt.ylabel("E(p)")
    # }}}
    # {{{ plot T1p data
    R1p = 1 / T1p
    fl.next(r"$R_{1}(p)$")
    fl.plot(R1p, "o", label="Experimental Data")
    # }}}
    # {{{Fit R1p
    T100 = 3.71
    dT10 = 0.2
    a, b, c, power = symbols("a b c power", real=True)  # symbols
    f = lmfitdata(R1p)  # initiate lmfit
    f.functional_form = (1 / (T100 + dT10 * power)) + (
        SL_conc_M / (a + b * power + c * power**2)
    )  # declare fitting function
    f.set_guess(
        a=dict(value=1.12e-4, min=1e-5, max=5e-3),
        b=dict(value=7.5e-5, min=-0.001, max=0.001),
        c=dict(value=-1.5e-5, min=-0.001, max=0.001),
    )
    f.settoguess()
    f.fit()
    R1p_fit = f.eval(Ep_pts)
    fl.plot(R1p_fit, ls=":", color="k", label="Fit", alpha=0.5)
    plt.ylabel(r"$R_{1} / s^{-1}$")
    plt.xlabel("Power / W")
    # }}}
    # {{{ Fit E(p)
    fl.next("E(p)")
    M0, A, phalf, power = symbols("M0 A phalf power", real=True)
    sp = power / (power + phalf)
    R1p = (1 / (T100 + dT10 * power)) + (
        SL_conc_M / (f.output("a") + f.output("b") * power + f.output("c") * power**2)
    )
    Ep_fit = lmfitdata(Ep["power", :-3])
    Ep_fit.functional_form = M0 - ((M0 * A * sp) / R1p)
    Ep_fit.set_guess(
        M0=dict(value=8e4, min=6.4e4, max=11e4),
        A=dict(value=3, min=0.5, max=17),
        phalf=dict(value=0.2, min=0.1, max=0.4),
    )
    Ep_fit.settoguess()
    Ep_fit.fit()
    thisfit = Ep_fit.eval(100)
    fl.plot(thisfit, ls=":", color="k", label="Fit", alpha=0.5)
    Ep_fit_text = r"M0 - $\frac{M0*A*s(p)}{R_{1}(p)}$"
    ksig = (Ep_fit.output("A") * ppt) / SL_conc_M
    ax = plt.gca()
    text(
        0.5,
        0.7,
        (3 * "\n") + list_symbs(Ep_fit),
        ha="center",
        va="center",
        size=10,
        transform=ax.transAxes,
    )
    text(
        0.5,
        0.75,
        r"$k_{\sigma}$ = %0.6f $M^{-1}s^{-1}$" % (ksig),
        ha="center",
        va="center",
        transform=ax.transAxes,
    )
    plt.xlabel("Power / W")
    plt.ylabel("E(p)")
    # }}}
