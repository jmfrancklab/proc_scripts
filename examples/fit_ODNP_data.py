"""Fitting ODNP Datasets for Ksigma
===================================
The T1(p) and E(p) integrals are stored in
previous post processing pulled from an H5 file.
The inverse of the T1(p) data processing is taken
to produce the R1(p) data which is fit using the
polyfit of the inverse of krho with two degrees of
freedom. The R1(p) fit is then fed into the
fitting routine for the enhancement data prior to
normalization. The cross relaxivity is then
calculated using the output of these fits as well
as the sample parameters that are fed to in (e.g.,
ppt value, and concentration).
"""
from pyspecdata import *
from sympy import symbols, Symbol, latex,lambdify
from scipy.io import loadmat

h5_file = "ras.h5"
nodename = "230706_M67_a"


# {{{ load data
Ep = find_file(f"{h5_file}", exp_type="AG_processed_data", expno=f"{nodename}/Ep")
# Some older h5 files save the T1p rather than the R1p. If there isn't an R1p expno then it will load the T1p integrals and convert to R1p by taking the inverse
try:
    R1p = find_file(f"{h5_file}", exp_type="AG_processed_data", expno=f"{nodename}/R1p")
except:
    T1p = find_file(f"{h5_file}", exp_type="AG_processed_data", expno=f"{nodename}/T1p")
    R1p = T1p**-1
# }}}
# {{{Find the index where the return powers begin
flip_idx = np.where(np.diff(Ep.fromaxis("power").data) < 0)[0][0] + 1
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
    krho_inv_coeff = krho_inv.polyfit("power", order=1)
    krho_inv_fine = R1p.fromaxis("power").eval_poly(krho_inv_coeff, "power")
    M0, A, phalf, p = symbols("M0 A phalf power", real=True)
    R1p_expression = (T10_p[0] + T10_p[1] * p) ** -1 + (
        Ep.get_prop("acq_params")["concentration"]
        / (krho_inv_coeff[0] + krho_inv_coeff[1] * p)
    )
    R1p_fit = lambdify(p,R1p_expression)
    fl.plot(R1p_fit(R1p.fromaxis('power')), ls=":", color="k", label="Fit", alpha=0.5)
    plt.ylabel(r"$R_{1} / s^{-1}$")
    plt.xlabel("Power / W")
    # }}}
    # {{{ Fit E(p)
    fl.next("Integrated Enhancement")
    sp = p / (p + phalf)
    Ep_fit = lmfitdata(Ep["power", :flip_idx])
    # Symbolic expression for Ep that is used in the symbolic function for the fitting of E(p)
    Ep_fit.functional_form = M0 - ((M0 * A * sp) / R1p_expression)
    A_guess = (
        1
        - (
            R1p["power", 0].data * (Ep["power", flip_idx].data / Ep["power", 0].data)
        ).real
    )
    Ep_fit.set_guess(
        M0=dict(value=Ep["power", 0].real.item(), min=1e4, max=11e4),
        A=dict(value=A_guess, min=0.5 * A_guess, max=3 * A_guess),
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
        (3 * "\n")
        + "\n".join(
            [f"${latex(Symbol(j))} = {k:0.5g}$" for j, k in Ep_fit.output().items()]
        ),
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
