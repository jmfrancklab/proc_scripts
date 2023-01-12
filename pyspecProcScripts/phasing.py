"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
import time
from matplotlib.patches import Ellipse
from matplotlib.transforms import blended_transform_factory
from scipy.optimize import minimize
from pylab import xlim, subplots, axvline, ylim, sca
import numpy as np
from numpy import r_, c_
from scipy import linalg
import logging
import matplotlib.pyplot as plt


def zeroth_order_ph(d, fl=None):
    r"""determine the covariance of the datapoints
    in complex plane, and use to phase the
    zeroth-order even if the data is both negative
    and positive.

    Parameters
    ==========
    d: nddata
        Complex data whose zeroth order phase you want
        to find.
    fl: figlist or None (default)
        If you want the diagnostic plots (showing the
        distribution of the data in the complex plane),
        set this to your figlist object.
        It will add a plot called "check covariance
        test"

    Returns
    =======
    retval: complex128
        The zeroth order phase of the data, as a
        complex number with magnitude 1.
        To correct the zeroth order phase of the data,
        divide by ``retval``.
    """
    cov_mat = np.cov(
        c_[d.data.real.ravel(), d.data.imag.ravel()].T,
        aweights=abs(d.data).ravel() ** 2  # when running proc_square_refl, having
        # this line dramatically reduces the size of the imaginary component
        # during the "blips," and the magnitude squared seems to perform
        # slightly better than the magnitude -- this should be both a robust
        # test and a resolution for issue #23
    )
    eigenValues, eigenVectors = linalg.eigh(cov_mat)
    mean_point = d.data.ravel().mean()
    mean_vec = r_[mean_point.real, mean_point.imag]
    # next 3 lines from stackexchange -- sort by
    # eigenvalue
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]  # first dimension x,y second evec #
    # determine the phase angle from direction of the
    # largest principle axis plus the mean
    # the vector in the direction of the largest
    # principle axis would have a norm equal to the
    # sqrt (std not variance) of the eigenvalue, except
    # that we only want to rotate when the distribution
    # is assymetric, so include only the excess of the
    # larger eval over the smaller
    assymetry_mag = float(sqrt(eigenValues[0]) - sqrt(eigenValues[1]))
    try:
        assym_ineq = (assymetry_mag * eigenVectors[:, 0] * mean_vec).sum()
    except:
        raise ValueError(
            strm(
                "check the sizes of the following:",
                size(assymetry_mag),
                size(eigenVectors),
                size(mean_vec),
            )
        )
    if assym_ineq > 0:
        # we want the eigenvector on the far side of the ellipse
        rotation_vector = mean_vec + assymetry_mag * eigenVectors[:, 0]
    else:
        rotation_vector = mean_vec - assymetry_mag * eigenVectors[:, 0]
    ph0 = np.arctan2(rotation_vector[1], rotation_vector[0])
    if fl is not None:
        d_forplot = d.C
        fl.next("check covariance test")
        fl.plot(
            d_forplot.data.ravel().real,
            d_forplot.data.ravel().imag,
            ".",
            alpha=0.25,
            label="before",
        )
        d_forplot /= np.exp(1j * ph0)
        fl.plot(
            d_forplot.data.ravel().real,
            d_forplot.data.ravel().imag,
            ".",
            alpha=0.25,
            label="after",
        )
        fl.plot(0, 0, "ko", alpha=0.5)
        fl.plot(mean_vec[0], mean_vec[1], "kx", label="mean", alpha=0.5)
        evec_forplot = (
            sqrt(eigenValues.reshape(1, 2)) * np.ones((2, 1)) * eigenVectors
        )  # scale by the std, not the variance!
        evec_forplot += mean_vec.reshape((-1, 1)) * np.ones((1, 2))
        fl.plot(
            evec_forplot[0, 0], evec_forplot[1, 0], "o", alpha=0.5, label="first evec"
        )
        fl.plot(evec_forplot[0, 1], evec_forplot[1, 1], "o", alpha=0.5)
        fl.plot(
            rotation_vector[1],
            "o",
            alpha=0.5,
            label="rotation vector",
        )
        norms = sqrt((evec_forplot ** 2).sum(axis=0))
        ell = Ellipse(
            xy=mean_vec,
            width=2 * sqrt(eigenValues[0]),
            height=2 * sqrt(eigenValues[1]),
            angle=180 / pi * np.arctan2(eigenVectors[1, 0], eigenVectors[0, 0]),
            color="k",
            fill=False,
        )
        ax = plt.gca()
        ax.set_aspect("equal", adjustable="box")
        ax.add_patch(ell)
    return np.exp(1j * ph0)


def ph1_real_Abs(s, dw, ph1_sel=0, ph2_sel=1, fl=None):
    r"""Performs first order phase correction with cost function
    by taking the sum of the absolute value of the real [DeBrouwer2009].

    .. todo::
        update with `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html>`_.

    Parameters
    ==========
    s: nddata
        Complex data whose first order phase you want
        to find.

    Returns
    =========
    retval

    .. rubric: References

    ..  [DeBrouwer2009] de Brouwer, H. (2009). Evaluation of algorithms for
        automated phase correction of NMR spectra. Journal of Magnetic
        Resonance (San Diego, Calif. : 1997), 201(2), 230–238.
        https://doi.org/10.1016/j.jmr.2009.09.017
    """
    if fl:
        fl.push_marker()
    ph1 = nddata(r_[-5:5:70j] * dw, "phcorr")
    dx = np.diff(ph1.getaxis("phcorr")[r_[0, 1]]).item()
    ph1 = np.exp(-1j * 2 * pi * ph1 * s.fromaxis("t2"))
    s_cost = s * ph1
    ph0 = s_cost.C.sum("t2")
    ph0 /= abs(ph0)
    s_cost /= ph0
    if fl:
        fl.next("phasing cost function")
    s_cost.run(np.real).run(abs).sum("t2")
    s_cost = s_cost["ph2", 1]["ph1", 0]
    if fl:
        fl.plot(s_cost, ".")
    s_cost.smoosh(["vd", "phcorr"], "transients")
    ph1_opt = np.asarray(s_cost.argmin("transients").item())
    logging.debug(strm("THIS IS PH1_OPT", ph1_opt))
    logging.debug(strm("optimal phase correction", repr(ph1_opt)))
    # }}}
    # {{{ apply the phase corrections
    def applyphase(arg, ph1):
        arg *= np.exp(-1j * 2 * pi * ph1 * arg.fromaxis("t2"))
        ph0 = arg.C.sum("t2")
        ph0 /= abs(ph0)
        arg /= ph0
        return arg

    def costfun(ph1):
        logging.debug(strm("PH1 IS THIS:", ph1))
        if type(ph1) is np.ndarray:
            ph1 = ph1.item()
        temp = s.C
        retval = applyphase(temp, ph1).run(np.real).run(abs).sum("t2").item()
        return retval

    r = minimize(costfun, ph1_opt, bounds=[(ph1_opt - dx, ph1_opt + dx)])
    assert r.success
    s = applyphase(s, r.x.item())
    fl.plot(r.x, r.fun, "x")
    fl.pop_marker()
    return s
    # }}}


def hermitian_function_test(
    s,
    direct="t2",
    aliasing_slop=3,  # will become a kwarg
    fl=None,
    basename=None,
    show_extended=False,
    echo_before=None,
):
    r"""determine the center of the echo via hermitian symmetry of the time domain.

    Parameters
    ==========
    direct:             str
                        Axis of data (i.e., direct dimension).
    aliasing_slop:      float
                        The multiple of the dwell time of the original 
                        signal that is sliced out due to the fact this amount of signal
                        is aliased at the beginning and end of the time-domain signal.
    basename:           str
                        Basename for plotting purposes of diagnostics
    show_extended:      bool
                        If true will show diagnostic of the zero filled signal
    echo_before:        float
                        1.5 times the echo time pulled from the acquisition parameters
    
    Returns
    =======
    echo_peak:          float
                        Time where the minimum of the cost function matches the peak of the echo
    cost_func_return    nddata
                        The cost function of the echo
    """
    # {{{ zero fill
    assert s.get_units(direct) is not None
    if s.get_ft_prop(direct):
        s_ext = s.C.ift(ift)
    else:
        s_ext = s.C
    t_dw = s.get_ft_prop(direct,'dt')
    print("trying to extend from",s_ext.getaxis(direct)[r_[0,-1]],"to",
            s_ext.getaxis(direct)[0] + 2 * np.diff(s_ext.getaxis(direct)[r_[0, -1]]).item())
    orig_bounds = s.getaxis(direct)[r_[0,-1]]
    orig_bounds -= orig_bounds[0] #starting point is at 0 for orig_bounds
    s_ext.extend(
        direct, s_ext.getaxis(direct)[0] + 2 * np.diff(s_ext.getaxis(direct)[r_[0, -1]]).item()
    )  # zero fill
    s_ext.register_axis({direct: 0}).ft(direct).set_ft_prop(
        direct, "start_time", 0
    ).ift(
        direct
    )  # forces the axis to *start* at 0
    if fl is not None:
        fl.push_marker()
        if basename is None:
            basename = f'randombasename{int(time.time()*10):d}'
        fl.basename = basename
        if show_extended:
            fl.next("data extended")
            if len(s_ext.dimlabels) > 1:
                fl.image(s_ext)
            else:
                fl.plot(s_ext, human_units=False)
        s_ext /= (
            abs(s_ext).mean_all_but(direct).data.max()
            )  # normalize by the average echo peak (for plotting purposes)
        fl.next("power terms")
        fl.plot(abs(s_ext).mean_all_but(
            direct)[direct:orig_bounds], label="echo envelope",
            human_units=False)
        fl.next("zoomed power terms")
        fl.plot(abs(s_ext).mean_all_but(
            direct)[direct:(0.1*orig_bounds)], label = "echo envelope",
            human_units = False)
    # }}}
    # {{{ the integral of the signal power up to t=Δt
    #     (first term in the paper)
    s_energy = s_ext.C
    s_energy.run(lambda x: abs(x) ** 2)
    s_energy.integrate(direct, cumulative=True)
    t_dwos = s_energy.get_ft_prop(direct, "dt")
    normalization_term = 2 * t_dwos / (s_energy.fromaxis(direct) + t_dwos)
    s_energy *= normalization_term
    s_energy.mean_all_but(direct)
    forplot = s_energy / t_dwos
    forplot.setaxis(direct, lambda x: x / 2)
    if fl is not None:
        fl.next("power terms")
        fl.plot(forplot[direct:orig_bounds], label="first energy term",
                human_units=False)
        fl.next("zoomed power terms")
        fl.plot(forplot[direct:(0.1*orig_bounds)], label="first energy term",
                human_units=False)
    # }}}
    # {{{ calculation the correlation between the echo and its hermitian
    #     conjugate
    #     → by definition, at the center of the echo, this should be
    #     equal to the previous (integral of the power) term
    s_correl = s_ext.C
    s_correl.ft(direct)
    s_correl.run(lambda x: x ** 2)
    s_correl.ift(direct)
    s_correl.mean_all_but(direct).run(abs)
    s_correl *= normalization_term
    forplot = s_correl / t_dwos
    forplot.setaxis(direct, lambda x: x / 2)
    if fl is not None:
        fl.next("power terms")
        fl.plot(forplot[direct:orig_bounds], label="correlation function",
                human_units=False)
        fl.next("zoomed power terms")
        fl.plot(forplot[direct:(0.1*orig_bounds)], label="correlation function",
                human_units=False)
    # }}}
    # {{{ calculate the cost function and determine where the center of the echo is!
    cost_func = s_energy - s_correl
    min_echo = aliasing_slop * t_dw
    cost_func.run(sqrt)  # based on what we'd seen previously (empirically), I
    #                     take the square root for a well-defined minimum -- it
    #                     could be better to do this before averaging in the
    #                     future
    forplot = cost_func / sqrt(t_dwos)
    forplot.setaxis(direct, lambda x: x / 2)
    cost_min = cost_func[direct:orig_bounds][direct:(min_echo,echo_before*2)].C.argmin(direct).item()
    if fl is not None:
        fl.next("power terms")
        fl.plot(
            forplot[direct:orig_bounds],
            label="cost function",
            c="violet",
            alpha=0.5,
            human_units=False,
        )
        fl.next("zoomed power terms")
        fl.plot(
            forplot[direct:(0.1*orig_bounds)],
            label="cost function",
            c="violet",
            alpha=0.5,
            human_units=False,
        )
    cost_func_return = forplot
    echo_peak = cost_min / 2.0
    if fl is not None:
        fl.next("power terms")
        fl.plot(echo_peak, forplot[direct:echo_peak].item(), "o", c="violet", alpha=0.3,
                human_units=False)
        axvline(x=echo_peak, linestyle=":")
        fl.next("zoomed power terms")
        fl.plot(echo_peak, forplot[direct:echo_peak].item(), "o", c="violet", alpha=0.3,
                human_units=False)
        axvline(x=echo_peak, linestyle=":")
        fl.pop_marker()
    # }}}
    return echo_peak,cost_func_return

def determine_sign(s, direct="t2", fl=None):
    """Given that the signal resides in `pathway`, determine the sign of the signal.
    The sign can be used, e.g. so that all data in an inversion-recover or
    enhancement curve can be aligned together.

    Parameters
    ==========
    s: nddata
        data with a single (dominant) peak, where you want to return the sign
        of the integral over all the data.
        This should only contain **a single coherence pathway**.
    direct: str (default "t2")
        Name of the direct dimension, along which the sum/integral is taken

    Returns
    =======
    data_sgn: nddata
        A dataset with all +1 or -1 (giving the sign of the original signal).
        Does *not* include the `direct` dimension
    """
    assert s.get_ft_prop(
        direct
    ), "this only works on data that has been FT'd along the direct dimension"
    if fl is not None:
        fl.push_marker()
        if basename is None:
            basename = f'randombasename{time.time()*10:d}'
        fl.basename = basename
        fl.next("selected pathway")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"))
    data_sgn = s.C.sum(direct)
    data_sgn /= zeroth_order_ph(data_sgn)
    data_sgn.run(np.real).run(lambda x: np.sign(x))
    if fl is not None:
        fl.next("check sign")
        fl.image((s.C.setaxis("vd", "#").set_units("vd", "scan #")) * data_sgn)
        fl.pop_marker()
    return data_sgn
