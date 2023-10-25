"""Fitting ODNP Datasets for Ksigma
===================================
The T1(p) and E(p) integrals are stored in previous post processing pulled from an H5 file. The inverse of the T1(p) data is taken to produce the R1(p) data which is fit using the polyfit of the inverse of krho with two degrees of freedom. The R1(p) fit is then fed into the fitting routine for the enhancement data prior to normalization. The cross relaxivity is then calculated using the output of these fits as well as the sample parameters that are fed to it (e.g., ppt value, and concentration).
"""
from pyspecdata import *
from sympy import symbols, Symbol, latex
from scipy.io import loadmat

target_directory = os.path.normpath(getDATADIR("AG_processed_data"))
h5_file = "ras.h5"
nodename = "230706_M67_a"
# {{{ plotting fn - for labels on final plot
def list_symbs(f):
    list_symbs = [f"${latex(Symbol(j))} = {k:0.5g}$" for j, k in f.output().items()]
    return "\n".join(list_symbs)


# }}}
# {{{ load data
Ep = nddata_hdf5(f"{h5_file}/{nodename}/Ep", directory=target_directory)
R1p = nddata_hdf5(f"{h5_file}/{nodename}/R1p", directory=target_directory)
# }}}
# {{{Find the index where the return powers begin
flip_idx = np.where(Ep.fromaxis("power").diff("power").data < 0)[0][0]
flip_idx += 1
# }}}
with figlist_var() as fl:
    # {{{Plot Ep
    fl.next("Integrated Enhancement")
    fl.plot(
        Ep["power", :flip_idx],
        "o",
        label="Progressive Saturation Data",
        capsize=6,
        alpha=0.5,
    )
    fl.plot(
        Ep["power", flip_idx:],
        "rx",
        label="Returning Power Check",
        capsize=6,
        alpha=0.5,
    )
    # }}}
    # {{{ plot R1p data
    fl.next(r"$R_{1}(p)$")
    fl.plot(R1p, "o", label="Experimental Data")
    # }}}
    # {{{load in T100 dataset
    T10_p = loadmat(
        search_filename(
            "T10_DI_water_230412", exp_type="AG_processed_data", unique=True
        )
    )["a"][0, :]
    R10_p = R1p.fromaxis("power").eval_poly(T10_p, "power") ** -1
    # }}}
    # {{{ fit krho inverse with two degrees of freedom and then apply to fit R1p
    krho_inv = Ep.get_prop("acq_params")["concentration"] / (R1p - R10_p)
    krho_inv_fit = krho_inv.polyfit("power", order=2)
    krho_inv_fine = R1p.fromaxis("power").eval_poly(krho_inv_fit, "power")
    R1p_fit = R10_p + Ep.get_prop("acq_params")["concentration"] / krho_inv_fine
    fl.plot(R1p_fit, ls=":", color="k", label="Fit", alpha=0.5)
    plt.ylabel(r"$R_{1} / s^{-1}$")
    plt.xlabel("Power / W")
    # }}}
    # {{{ Fit E(p)
    fl.next("Integrated Enhancement")
    M0, A, phalf, p = symbols("M0 A phalf power", real=True)
    sp = p / (p + phalf)
    R1p = (T10_p[0] + T10_p[1] * p) ** -1 + (
        Ep.get_prop("acq_params")["concentration"]
        / (krho_inv_fit[0] + krho_inv_fit[1] * p + krho_inv_fit[2] * p**2)
    )  # Symbolic expression for R1p that is used in the symbolic function for the fitting of E(p)
    Ep_fit = lmfitdata(Ep["power", :flip_idx])
    Ep_fit.functional_form = M0 - ((M0 * A * sp) / R1p)
    Ep_fit.set_guess(
        M0=dict(value=Ep["power", 0].real.item(), min=1e4, max=11e4),
        A=dict(value=13, min=0.5, max=17),
        phalf=dict(value=0.2, min=0.1, max=0.4),
    )
    Ep_fit.settoguess()
    Ep_fit.fit()
    thisfit = Ep_fit.eval(100)
    fl.plot(thisfit, ls=":", color="k", label="Fit", alpha=0.5)
    ksig = (
        Ep_fit.output("A") * Ep.get_prop("acq_params")["guessed_MHz_to_GHz"] * 1e-3
    ) / Ep.get_prop("acq_params")["concentration"]
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
    plt.ylabel(r"$M_{0}E(p)$")
    # }}}
