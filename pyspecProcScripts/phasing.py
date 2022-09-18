"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
import time
from matplotlib.patches import Ellipse
from matplotlib.transforms import blended_transform_factory
from scipy.optimize import minimize
from pylab import xlim, subplots, axvline, ylim, sca, rand
import numpy as np
from numpy import r_, c_
from scipy import linalg
import scipy.signal.windows as sci_win
import logging
import matplotlib.pyplot as plt

def det_devisor(fl):
    if fl.units[fl.current] == "s":
        divisor = 1
    elif fl.units[fl.current] == "ms":
        divisor = 1e-3
    elif fl.units[fl.current] == "\\mu s":
        divisor = 1e-6
    else:
        raise ValueError(
            f"current units are {fl.units[fl.current]} right now, only programmed to work with results of s, ms and μs"
        )
    return divisor

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
            rotation_vector[1], "o", alpha=0.5, label="rotation vector",
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
    amp_threshold=0.05,  # region over which we have signal
    fl=None,
    basename=None,
    show_extended=False,
    upsampling=40,
    energy_threshold=0.98,
    energy_threshold_lower=0.1,
    enable_refinement=False,
):
    r"""Determine the center of the echo via hermitian symmetry of the time domain.

    Note the following issues/feature:

    -   Data *must* start at t=0 (due to the way that we calculate the energy
        term in the expression)
    -   Data typically has a rising edge, which can lead to ringing.
        Therefore, we apply a Tukey filter to the data.  Because this
        suppresses the edges in frequency space, you will likely want to slice
        a slightly wider bandwidth than you are interested in.
    -   If you have an initial portion of signal that has very low SNR (i.e. if
        your echo time is long relative to :math:`T_2^*`), that will affect the
        SNR of the cost function.  Trailing noise will not affect the SNR of
        the cost function, only lead to slightly longer calculation times.

    Parameters
    ==========
    direct:             str
        Axis of data (i.e., direct dimension).
    energy_threshold:   float
        To avoid picking a minimum that's based off of very little data, only
        recognize the cost function over an interval where the normalized
        energy function is equal to `energy_threshold` of its max value.
    aliasing_slop: int (default 3)
        The signal typically has an initial rising
        portion where -- even due to FIR or
        instrumental filters in the raw data,
        the signal has some kind of Gibbs ringing in the
        transition between signal an zero.
        This portion of the signal will not mirror the
        signal on the other side of the echo, and
        should not even be included in the cost
        function!
        We believe the most consistent way to specify
        the amount of signal that should be excluded as
        an number of datapoints, hence why we ask for
        an integer here.

    .. todo::

        AG fix docstring + cite coherence paper
    """
    # {{{ zero fill and upsample
    # {{{ get in time domain and grab dwell time
    assert s.get_units(direct) is not None
    if s.get_ft_prop(direct):
        s_timedom = s.C.ift(direct)
    else:
        s_timedom = s.C
    s_ext = s_timedom.C
    assert (
        s.getaxis(direct)[0] == 0.0
    ), """In order to
    calculate the signal energy term correctly, the
    signal must start at t=0  so set the start of the
    acquisition in the *non-aliased* time domain to 0 (something like
    data['t2'] -= acqstart) to avoid confusion"""
    t_dw = s.get_ft_prop(direct, "dt")
    orig_bounds = s.getaxis(direct)[r_[0, -1]]
    plot_bounds = orig_bounds # allow this to be set to ± inf if I want to see all
    # }}}
    # {{{ force the axis to *start* at 0
    #     since we're doing FT here, also extend to
    #     twice the length!
    logging.debug(strm("initial size", ndshape(s_ext), "direct is", direct))
    s_ext.ft(
        direct, pad=2 ** int(np.ceil(np.log(ndshape(s_ext)[direct]) / np.log(2)) + 1)
    )
    assert (
        s_ext.get_ft_prop(direct, "start_time") == 0
    ), "FT start point should also be equal to zero -- doing otherwise doesn't make sense"
    # }}}
    # {{{ move into over-sampled time domain -- do this
    # *after* previous to avoid aliasing glitch
    tukeyfilter = s_ext.fromaxis(direct).run(lambda x: sci_win.tukey(len(x)))
    s_ext *= tukeyfilter
    s_ext.ift(direct, pad=2 ** int(np.ceil(np.log(ndshape(s)[direct] * upsampling) / np.log(2))))
    s_ext[direct:(orig_bounds[-1],None)] = 0 # explicitly zero, in case there are aliased negative times!
    # }}}
    # {{{ now I need to throw out the initial, aliased
    #     portion of the signal -- do this manually by
    #     index -- this is more lines than needed in
    #     order to be explanatory.  Note that I
    #     plot the signal before I actually throw stuff
    #     out
    t_dwos = s_ext.get_ft_prop(direct, "dt") # oversampled dwell
    min_echo = aliasing_slop * t_dw
    min_echo_idx = int(min_echo/t_dwos + 0.5)
    min_echo = min_echo_idx * t_dwos
    if fl is not None:
        fl.push_marker()
        if basename is None:
            basename = f"randombasename{int(rand()*1e5):d}"
        fl.basename = basename
        if show_extended:
            fl.next("data extended")
            if len(s_ext.dimlabels) > 1:
                fl.image(s_ext)
            else:
                fl.plot(s_ext)
        s_ext /= (
            abs(s_ext).mean_all_but(direct).data.max()
        )  # normalize by the average echo peak (for plotting purposes)
        fl.next("power terms")
        forplot = abs(s_ext).mean_all_but(direct)[direct:plot_bounds]
        forplot[direct] -= min_echo # so zero is the
        #                             first part of the
        #                             echo after the
        #                             aliasing slop, as
        #                             indicated below
        fl.plot(
            forplot, label="echo envelope",
        )
    s_ext[direct, :-min_echo_idx] = s_ext[direct, min_echo_idx:]
    if fl is not None:
        fl.next("power terms")
        forplot = abs(s_ext).mean_all_but(direct)[direct:plot_bounds]
        fl.plot(
            forplot, label="echo envelope",
        )
    # }}}
    # }}}
    # {{{ the integral of the signal power up to t=Δt
    #     (first term in the paper)
    s_energy = s_ext.C
    s_energy.run(lambda x: abs(x) ** 2)
    s_energy.integrate(direct, cumulative=True)
    s_energy.mean_all_but(direct)
    t_dwos = s_energy.get_ft_prop(direct, "dt")
    normalization_term = 2 * t_dwos / (s_energy.fromaxis(direct) + t_dwos)
    s_energy *= normalization_term
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
    # }}}
    # {{{ calculate the cost function and determine where the center of the echo is!
    cost_func = abs(s_energy - s_correl) # b/c this should not be less than 0, so penalize for numerical error when it's not!
    reasonable_energy_range = s_energy.contiguous(lambda x: abs(x) > energy_threshold*abs(x.data).max())[0,:]
    _,reasonable_energy_range[1] = s_energy.contiguous(lambda x: abs(x) >
            energy_threshold*energy_threshold_lower*abs(x.data).max())[0,:]
    print(reasonable_energy_range)
    cost_func = cost_func[direct:reasonable_energy_range]
    cost_func.run(lambda x: x/sqrt(abs(x)))  # based on what we'd seen
    #             previously (empirically), I take the
    #             square root for a well-defined
    #             minimum -- it could be better to do
    #             this before averaging in the future
    cost_min = cost_func.C.argmin(direct).item()
    if fl is not None:
        forplot = s_energy / t_dwos
        forplot.mean_all_but(direct)
        forplot.setaxis(direct, lambda x: x / 2)
        fl.plot(forplot[direct:plot_bounds], label="first energy term")
        forplot = s_correl / t_dwos
        forplot.mean_all_but(direct)
        forplot.setaxis(direct, lambda x: x / 2)
        fl.plot(forplot[direct:plot_bounds], label="correlation function")
        forplot = cost_func / sqrt(t_dwos)
        forplot.setaxis(direct, lambda x: x / 2)
        fl.plot(
            forplot,
            label="cost function",
            c="violet",
            alpha=0.5,
        )
    echo_peak = cost_min / 2.0
    if fl is not None:
        fl.plot(
            echo_peak / det_devisor(fl),
            forplot[direct:echo_peak].item(),
            "o",
            c="violet",
            alpha=0.3,
            human_units=False,
        )
        axvline(x=echo_peak / det_devisor(fl), linestyle=":")
    # }}}
    echo_idx = int((echo_peak + min_echo)/t_dw+0.5)
    if enable_refinement and echo_idx-aliasing_slop > 3:
        s_foropt = s_timedom[direct,
                aliasing_slop:2*(echo_idx-aliasing_slop)+1]
        center_idx = (2*(echo_idx-aliasing_slop)+1-aliasing_slop)//2+1
        if fl is not None:
            fl.next('foropt')
            fl.plot(abs(s_foropt), human_units=False)
            axvline(x=(echo_peak + min_echo), alpha=0.5)
            axvline(x=s_foropt.getaxis(direct)[ndshape(s_foropt)[direct]//2+1],
                    color='r',ls=':', alpha=0.5)
        shifts = nddata(
                t_dw*r_[-1.5:1.5:300j],'echo shift')
        s_foropt.ft(direct)
        # negative pushes time domain to right --
        # should represent that echo occured earlier
        # than expected, needed to be pushed right to
        # be centered
        s_foropt *= np.exp(-1j*2*np.pi*shifts*s_foropt.fromaxis(direct))
        s_foropt.ift(direct)
        ph0 = s_foropt[direct,
                center_idx].C
        ph0 = zeroth_order_ph(ph0)
        s_foropt /= ph0
        s_foropt = s_foropt[direct,2:-2]
        s_foropt -= s_foropt.C[direct,::-1].run(np.conj)
        s_foropt.run(lambda x: abs(x)**2)
        s_foropt.mean_all_but('echo shift')
        if fl is not None:
            fl.next('refinement')
            fl.plot(s_foropt, human_units=False)
        s_foropt = s_foropt.argmin('echo shift').item()
        if fl is not None:
            fl.next('refinement')
            axvline(x=s_foropt)
            fl.next("power terms")
            axvline(x=(t_dw*echo_idx-s_foropt -min_echo)/det_devisor(fl),
                    ls='--',
                    alpha=0.25)
        if fl is not None:
            fl.pop_marker()
        return t_dw*echo_idx - s_foropt
    else:
        logging.info("warning: can't do hermitian phasing refinement -- not enough points")
        if fl is not None:
            fl.pop_marker()
        return echo_peak + min_echo
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
            basename = f"randombasename{int(rand()*1e5):d}"
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
