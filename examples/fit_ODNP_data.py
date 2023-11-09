"""Fitting ODNP Datasets for Ksigma
=================================== 
The T1(p) and E(p) integrals are generated
using the example, `generate_integrals.py`
(as of 11/8/23 the `generate_integrals` is not
updated), and stored in
previous post processing pulled from an H5 file.
The inverse of the :math:`T_1(p)` data processing is taken
to produce the :math:`R_1(p)` data which is fit using the
polyfit of the inverse of krho with two degrees of
freedom. The :math:`R_1(p)` fit is then fed into the
fitting routine for the enhancement data prior to
normalization. The cross relaxivity is then
calculated using the output of these fits as well
as the sample parameters that are fed to in (*e.g.*,
ppt value, and concentration).
"""
from pyspecdata import *
from sympy import symbols, Symbol, latex, lambdify
from scipy.io import loadmat

# {{{ This block changes when the data changes -- everything else should be
# left alone for most circumstances
data_info = dict(
    filename="ras.h5",  # h5 file containing table of integrals for different datasets
    data_dir="AG_processed_data",  # directory of the dataset of table of integrals
    nodename="230706_M67_a",
)  # specific nodename of the dataset of interest
T100_info = dict(
    filename="T10_DI_water_230412",
    data_dir="AG_processed_data",
)
# }}}

# {{{ load data
Ep = find_file(
    data_info["filename"],
    exp_type=data_info["data_dir"],
    expno=f"{data_info['nodename']}/Ep",
)
# {{{ Some older h5 files save the T1p rather than the R1p. If there isn't an
# R1p expno then it will load the T1p integrals and convert to R1p by taking
# the inverse
try:
    R1p = find_file(
        data_info["filename"],
        exp_type=data_info["data_dir"],
        expno=f"{data_info['nodename']}/R1p",
    )
except:
    T1p = find_file(
        data_info["filename"],
        exp_type=data_info["data_dir"],
        expno=f"{data_info['nodename']}/T1p",
    )
    R1p = 1 / T1p
# }}}
# }}}
# {{{ The powers go up, and then go back down in order to check for
# reproducibility.  Figure out where this flip
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
    # {{{ fit kᵨ⁻¹ with two degrees of freedom (to a straight line) and then
    # apply to fit R1p
    T10_p = loadmat(
        search_filename(
            T100_info["filename"], exp_type=T100_info["data_dir"], unique=True
        )
    )["a"][0, :]
    R10_p = 1 / (R1p.fromaxis("power").eval_poly(T10_p, "power"))
    powers_fine = nddata(
        r_[0 : R1p.getaxis("power") : 300j], "p"
    )  # PR COMMENT: related to my comments below
    krho_inv = Ep.get_prop("acq_params")["concentration"] / (R1p - R10_p)
    krho_inv_coeff = krho_inv.polyfit("power", order=1)
    krho_inv_fine = R1p.fromaxis("power").eval_poly(
        krho_inv_coeff, "power"
    )  # PR COMMENT: why is this called "fine" it's only evaluated at the powers present in the R1p data, which is not that many
    M0, A, phalf, p = symbols("M0 A phalf power", real=True)
    R1p_expression = (T10_p[0] + T10_p[1] * p) ** -1 + (
        Ep.get_prop("acq_params")["concentration"]
        / (krho_inv_coeff[0] + krho_inv_coeff[1] * p)
    )
    R1p_fit = lambdify(p, R1p_expression)
    fl.plot(
        R1p_fit(R1p.fromaxis("power")),
        ls=":",  # PR COMMENT: again, why are you plotting this as :, since it's only evaluated at the powers of interest -- maybe you want to construct a new axis?  Why not just plot as a solid line, anyways?
        color="k",
        label="Fit",
        alpha=0.5,
    )
    plt.ylabel(r"$R_{1} / s^{-1}$")
    plt.xlabel("Power / W")
    # }}}
    # {{{ Fit E(p)
    fl.next("Integrated Enhancement")
    sp = p / (p + phalf)
    Ep_fit = lmfitdata(Ep["power", :flip_idx])
    # Symbolic expression for Ep that is used in the symbolic function for the fitting of E(p)
    Ep_fit.functional_form = M0 - ((M0 * A * sp) / R1p_expression)
    # generate a guess for the A parameter of the fit based on the normalized
    # enhancement weighted by the relaxation rate. The bounds for the fit are
    # then set to center around this value.
    # PR COMMENT: what's wrong with this reasoning?
    # Since
    # E(pₘₐₓ) = 1 - A s(pₘₐₓ)/R₁(pₘₐₓ)
    # and the max s(p) is about 1,
    # 1-E(pₘₐₓ)R₁(pₘₐₓ) ≈ A
    A_guess = (
        1
        - (
            R1p["power", 0].real.item()
            * (
                Ep["power", flip_idx].real.item()
            )  # PR COMMENT: you were dividing by 1?
        ).real
    )
    Ep_fit.set_guess(
        M0=dict(value=Ep["power", 0].real.item(), min=1e4, max=11e4),
        A=dict(value=A_guess, min=0.2 * A_guess, max=3 * A_guess),
        phalf=dict(value=0.2, min=0.05, max=1.0),
    )
    Ep_fit.settoguess()
    Ep_fit.fit()
    thisfit = Ep_fit.eval(100)
    fl.plot(thisfit, ls=":", color="k", label="Fit", alpha=0.5)
    ksig = (
        Ep_fit.output("A")
        * Ep.get_prop("acq_params")["guessed_MHz_to_GHz"]
        * 1e-3  # guessed_MHz_to_GHz is the actual acquired ppt but the config
        # file currently does not have a key for actual ppt value
        # PR COMMENT: from what you were telling me, it just overwrites this value? if so, say so here.
    ) / Ep.get_prop("acq_params")["concentration"]
    ax = plt.gca()
    text(
        0.5,
        0.7,
        (3 * "\n")
        + "\n".join(
            [
                f"${latex(Symbol(j))} = {k:0.5g}$"
                for j, k in Ep_fit.output().items()
            ]
            + [r"$k_{\sigma} = %0.6f M^{-1}s^{-1}$" % (ksig)]
        ),
        ha="center",
        va="center",
        size=10,
        transform=ax.transAxes,
    )
    plt.xlabel("Power / W")
    plt.ylabel(r"$M_{0}E(p)$")
    # }}}
