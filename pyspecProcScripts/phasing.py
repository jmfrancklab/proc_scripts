"""This module includes routines for phasing NMR spectra."""

from pyspecdata import nddata, ndshape, strm
from matplotlib.patches import Ellipse
from scipy.optimize import minimize
from pylab import axvline, legend
import numpy as np
from numpy import r_, sqrt, pi
from scipy import linalg
import scipy.signal.windows as sci_win
import logging
import matplotlib.pyplot as plt
from .simple_functions import select_pathway
from itertools import cycle

default_matplotlib_cycle = cycle(
    plt.rcParams["axes.prop_cycle"].by_key()["color"]
)


def det_devisor(fl):
    if fl.units[fl.current] == "s":
        divisor = 1
    elif fl.units[fl.current] == "ms":
        divisor = 1e-3
    elif fl.units[fl.current] == "\\mu s":
        divisor = 1e-6
    else:
        raise ValueError(
            f"current units are {fl.units[fl.current]} right now, only"
            " programmed to work with results of s, ms and μs"
        )
    return divisor


def zeroth_order_ph(d, fl=None):
    r"""determine the moment of inertial of the datapoints
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
    realvector = d.data.real.ravel()
    imvector = d.data.imag.ravel()
    R2 = np.mean(realvector**2)
    I2 = np.mean(imvector**2)
    C = np.mean(
        realvector * imvector
    )  # for moment of inertia, this term is negative, but we use positive
    # instead, so that the ellipse is aligned with the distribution
    # note that this effectively changes the relative sign of x and y,
    # reflecting the ellipse so that it circles the elements rather than going
    # around them, and it's note the same as just flipping the eigenvalues
    inertia_matrix = np.array(
        [[R2, C], [C, I2]]
    )  # moment of inertia, with C inverted -- see comment below
    eigenValues, eigenVectors = linalg.eigh(inertia_matrix)
    mean_point = d.data.ravel().mean()
    # next 3 lines from stackexchange -- sort by
    # eigenvalue
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]
    # eigenVectors[1,:] *= -1   # leave this line -- uncommenting this and
    #                             negating the C above yields the same result!
    rotation_vector = eigenVectors[:, 0]
    ph0 = np.arctan2(rotation_vector[1], rotation_vector[0])
    if fl is not None:
        fl.push_marker()
        d_forplot = d.C
        fl.next("check covariance test")
        fl.plot(
            d_forplot.data.ravel().real,
            d_forplot.data.ravel().imag,
            ".",
            alpha=0.25,
            label="before",
        )
        plt.xlabel("real")
        plt.ylabel("imag")
        d_forplot /= np.exp(1j * ph0)
        fl.plot(
            d_forplot.data.ravel().real,
            d_forplot.data.ravel().imag,
            ".",
            alpha=0.25,
            label="after",
        )
        fl.plot(0, 0, "ko", alpha=0.5)
        evec_forplot = (
            2
            * sqrt(eigenValues.reshape(1, 2))
            * np.ones((2, 1))
            * eigenVectors
        )  # scale by the std, not the variance!
        fl.plot(
            evec_forplot[0, 0],
            evec_forplot[1, 0],
            "o",
            alpha=0.5,
            label="first evec",
        )
        fl.plot(evec_forplot[0, 1], evec_forplot[1, 1], "o", alpha=0.5)
        ell = Ellipse(
            xy=[0, 0],
            width=4 * sqrt(eigenValues[0]),
            height=4 * sqrt(eigenValues[1]),
            angle=180
            / pi
            * np.arctan2(eigenVectors[1, 0], eigenVectors[0, 0]),
            color="k",
            fill=False,
        )
        ax = plt.gca()
        ax.set_aspect("equal", adjustable="box")
        ax.add_patch(ell)
        fl.pop_marker()
    return np.exp(1j * ph0) * np.sign((mean_point / np.exp(1j * ph0)).real)


def ph1_real_Abs(s, dw, ph1_sel=0, ph2_sel=1, fl=None):
    r"""Performs first order phase correction with cost function
    by taking the sum of the absolute value of the real [DeBrouwer2009].

    .. todo::
        update with `sphinxcontrib-bibtex
        <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html>`_.

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


def fid_from_echo(
    d,
    signal_pathway,
    fl=None,
    add_rising=False,
    direct="t2",
    exclude_rising=3,
    slice_multiplier=20,
    peak_lower_thresh=0.1,
    show_hermitian_sign_flipped=False,
    show_shifted_residuals=False,
    frq_center=None,
    frq_half=None,
):
    """

    Parameters
    ==========
    signal_pathway: dict
        coherence transfer pathway that correspond to the signal
    fl: figlist or None (default)
        If you want the diagnostic plots (showing the distribution of the
        data in the complex plane), set this to your figlist object.
    add_rising: boolean
        Take the first part of the echo (that which rises to the maximum)
        and add the decaying (FID like) part. This increases the SNR of
        the early points of the signal.
    direct: string
        Name of the direct dimension
    exclude_rising: int
        In general it is assumed that the first few points of signal
        might be messed up due to dead time or other issues (assuming a
        tau of 0).  This option allows us to add a rising edge to the
        echo and exclude the first few points of the signal. Note, to use
        this, add_rising must be True.
    slice_multiplier: int
        The calculated frequency slice is calculated by taking the center
        frequency and extending out to values that are included in the
        peak times this multiplier. Therefore the larger this value the
        larger the frequency slice. Increasing this value might serve
        useful in the case of noisy spectra.
    peak_lower_thresh: float
        Fraction of the signal intensity used in calculating the
        frequency slice. The smaller the value, the wider the slice.
    show_hermitian_sign_flipped: boolean
        Diagnostic in checking the sign of the signal prior to the
        hermitian phase correction
    show_shifted_residuals: boolean
        Diagnostic in analyzing the residuals after the hermitian phase
        correction.
    frq_center: float (default None)
        The center of the peak.
        This only exists so that we don't end up calling
        `find_peakrange` redundantly,
        and it should come from a previous call to `find_peakrange` if
        it's used.
    frq_half: float (default None)
        The half-width of the peak.
        This only exists so that we don't end up calling
        `find_peakrange` redundantly,
        and it should come from a previous call to `find_peakrange` if
        it's used.

    Returns
    =======
    d: nddata
        FID of properly sliced and phased signal
    """
    if frq_center is None:
        frq_center, frq_half = find_peakrange(
            d, fl=fl, direct=direct, peak_lower_thresh=peak_lower_thresh
        )
    if fl is not None and "autoslicing!" in fl:
        fl.next("autoslicing!")
        axvline(
            x=frq_center - slice_multiplier * frq_half,
            color="k",
            ls="--",
            alpha=0.5,
            label="final slice",
        )
        axvline(
            x=frq_center + slice_multiplier * frq_half,
            color="k",
            ls="--",
            alpha=0.5,
            label="final slice",
        )
        legend()
    slice_range = r_[-1, 1] * slice_multiplier * frq_half + frq_center
    reduced_slice_range = r_[-1, 1] * 2 * frq_half + frq_center
    # }}}
    d = d[direct:slice_range]
    d.ift(direct)
    # {{{ apply phasing, and check the residual
    d[direct] -= d.getaxis(direct)[0]
    if fl is not None and fl.basename is not None:
        thebasename = fl.basename
    else:
        thebasename = ""
    # {{{ sign flip and average input for hermitian
    input_for_hermitian = select_pathway(d, signal_pathway).C
    signflip = input_for_hermitian.C.ft(direct)[direct:reduced_slice_range]
    idx = abs(signflip).mean_all_but(direct).data.argmax()
    signflip = signflip[direct, idx]
    ph0 = zeroth_order_ph(signflip)
    signflip /= ph0
    signflip.run(np.real)
    signflip /= abs(signflip)
    input_for_hermitian /= signflip
    if fl is not None and show_hermitian_sign_flipped:
        fl.next("sign flipped for hermitian")
        input_for_hermitian.reorder(direct, first=False)
        fl.image(input_for_hermitian)
    input_for_hermitian.mean_all_but(direct)
    # }}}
    best_shift = hermitian_function_test(
        input_for_hermitian,
        basename=" ".join([thebasename, "hermitian"]),
        fl=fl,
    )
    d_save = d.C
    t_dw = d_save.get_ft_prop(direct, "dt")
    if show_shifted_residuals:
        test_array = (
            t_dw / 3 * r_[-6, -3, -1, 1, 3, 6, 0]
        )  # do 0 last, so that's what it uses
    else:
        test_array = r_[0]
    for test_offset in test_array:
        test_shift = best_shift + test_offset
        d = d_save.C
        if test_offset == 0:
            thiscolor = "k"
            zeroth_fl = fl
        else:
            thiscolor = next(default_matplotlib_cycle)
            zeroth_fl = None
        if show_shifted_residuals:
            d.set_plot_color(thiscolor)
        d.setaxis(direct, lambda x: x - test_shift).register_axis({direct: 0})
        ph0 = zeroth_order_ph(
            select_pathway(d, signal_pathway)[direct:0.0], fl=zeroth_fl
        )
        d /= ph0
        if fl is not None:
            t_start = d.getaxis(direct)[0]
            d_sigcoh = select_pathway(d, signal_pathway)[
                direct : (t_start, -2 * t_start)
            ]
            d_sigcoh = select_pathway(d, signal_pathway).squeeze()
            s_flipped = d_sigcoh[direct:(t_start, -t_start)][direct, ::-1].C
            idx = (ndshape(s_flipped)[direct]) // 2
            ph0 = zeroth_order_ph(d_sigcoh[direct, idx])
            d_sigcoh /= ph0
            s_flipped /= ph0
            s_flipped.run(np.conj)
            s_flipped[direct] = s_flipped[direct][::-1]
            if (
                s_flipped.getaxis(direct)[idx] == 0
                and d_sigcoh.getaxis(direct)[idx] == 0
            ):  # should be centered about zero,
                # but will not be if too lopsided
                for_resid = (
                    abs(s_flipped - d_sigcoh[direct:(t_start, -t_start)]) ** 2
                )
                N_ratio = for_resid.data.size
                for_resid.mean_all_but(direct).run(sqrt)
                N_ratio /= for_resid.data.size  # the signal this
                #                                  has been plotted
                #                                  against is signal
                #                                  averaged by
                #                                  N_ratio
                fl.next("residual after shift")
                fl.plot(
                    for_resid / sqrt(N_ratio),
                    human_units=False,
                    label="best shift%+e, mean of residual" % test_offset,
                )
                # {{{ zoom to twice the width of the displayed residual
                #     on the left and more on the right
                left_bound = for_resid[direct][0]
                right_bound = for_resid[direct][-1]
                left_bound, right_bound = (
                    r_[-1.0, 4.0] * (right_bound - left_bound)
                    + (right_bound + left_bound) / 2
                )
                # }}}
                if test_offset == 0:
                    fl.plot(
                        d_sigcoh.C.mean_all_but(direct).run(abs),
                        alpha=0.8,
                        human_units=False,
                        label="best shift%+e, abs of mean" % test_offset,
                    )
                    s_flipped.set_plot_color("r")
                    fl.plot(
                        s_flipped.C.mean_all_but(direct).run(abs),
                        alpha=0.5,
                        human_units=False,
                        label="best shift%+e, abs of flipped mean"
                        % test_offset,
                    )
                    ax = plt.gca()
                    yl = ax.get_ylim()
                    ax.set_ylim((0, yl[-1]))
                plt.gca().set_xlim((left_bound, right_bound))
    # }}}
    if add_rising:
        d_rising = d[direct:(None, 0)][
            direct, exclude_rising:-1
        ]  # leave out first few points
        d_rising.run(np.conj)
        N_rising = ndshape(d_rising)[direct]
    d = d[direct:(0, None)]
    if add_rising:
        d[direct, 1 : N_rising + 1] = (
            d[direct, 1 : N_rising + 1] + d_rising[direct, ::-1]
        ) / 2
    d *= 2
    d[direct, 0] *= 0.5
    d.ft(direct)
    return d


def find_peakrange(
    d,
    direct="t2",
    peak_lower_thresh=0.1,
    peak_lowest_thresh=0.03,
    width_guess=50.0,
    max_echo=20e-3,
    fl=None,
):
    """find the range of frequencies over which the signal occurs, so that we
    can autoslice.
    Always assume that the signal is symmetric about zero and confined to the central 1/2 of the spectrum.

    Parameters
    ==========
    d : nddata
        Data in the frequency domain -- will not be altered.
    direct : str (default "t2")
        The name of the direct dimension
    peak_lower_thresh: float
        Fraction of the signal intensity used in calculating the
        frequency slice. The smaller the value, the wider the slice.
    width_guess: float
        Guess of the signal width, in Hz.
    max_echo: float
        Echo can be away from zero by this much (in s).
    fl : figlist (default None)
        If you want to see diagnostic plots, feed the figure list.

    Returns
    =======
    frq_center : float
        The midpoint of the frequency slice.
    frq_half : float
        Half the width of the frequency slice.
        Given in this way, so you can easily do
        >>> newslice = r_[-expansino,expansion]*frq_half+frq_center
    """
    # {{{ autodetermine slice range
    freq_envelope = d.C
    freq_envelope.ift(direct)
    freq_envelope = freq_envelope[
        direct:(0, None)
    ]  # slice out rising echo estimate according to experimental tau in order
    #   to limit oscillations
    freq_envelope[direct, 0] *= 0.5
    freq_envelope.ft(direct)
    freq_envelope.mean_all_but(direct).run(abs)
    if fl is not None:
        fl.next("autoslicing!")
        fl.plot(freq_envelope, human_units=False, label="signal energy")
    freq_envelope.convolve(
        direct,
        width_guess,
        enforce_causality=False,
    )
    SW = 1 / freq_envelope.get_ft_prop(direct, "dt")
    # baseline, assuming we're constrained to the middle 1/3 of the SW
    freq_envelope -= (
        freq_envelope[direct : tuple(-r_[0.165, 0.18] * SW)].mean().item()
        + freq_envelope[direct : tuple(r_[0.165, 0.18] * SW)].mean().item()
    ) / 2
    if fl is not None:
        fl.next("autoslicing!")
        fl.plot(
            freq_envelope,
            human_units=False,
            label="signal energy\nconv + baselined",
        )
    narrow_ranges = freq_envelope.contiguous(lambda x: x > 0.5 * x.data.max())
    wide_ranges = freq_envelope.contiguous(
        lambda x: x > peak_lower_thresh * x.data.max()
    )
    widest_ranges = freq_envelope.contiguous(
        lambda x: x > peak_lowest_thresh * x.data.max()
    )

    def filter_ranges(B, A):
        """where A and B are lists of ranges (given as tuple pairs), filter B
        to only return ranges that include ranges given in A"""
        return [
            np.array(b)
            for b in B
            if any(b[0] <= a[0] and b[1] >= a[1] for a in A)
        ]

    peakrange = filter_ranges(wide_ranges, narrow_ranges)
    peakrange = filter_ranges(widest_ranges, peakrange)
    if len(peakrange) > 1:
        max_range_width = max(
            [thisrange[1] - thisrange[0] for thisrange in peakrange]
        )
        range_gaps = [
            peakrange[j + 1][0] - peakrange[j][1]
            for j in range(len(peakrange) - 1)
        ]
        # {{{ if the gaps are all smaller than the max peak that was found, we
        #     just have "breaks" in the peak, so merge them.  Otherwise, fail.
        if any(np.array(range_gaps) > max_range_width):
            if fl is not None:
                fl.next("debug filter ranges")
                fl.plot(freq_envelope, human_units=False)
                for thisrange in peakrange:
                    fl.plot(freq_envelope[direct:thisrange], human_units=False)
            raise ValueError("finding more than one peak!")
        else:
            peakrange = [(peakrange[0][0], peakrange[-1][1])]
    assert len(peakrange) == 1
    peakrange = peakrange[0]
    # }}}
    frq_center = np.mean(peakrange).item()
    frq_half = np.diff(peakrange).item() / 2
    d.set_prop("peakrange", peakrange)
    if fl is not None:
        fl.next("autoslicing!")
        axvline(x=frq_center, color="k", alpha=0.5, label="center frq")
        axvline(
            x=frq_center - frq_half,
            color="k",
            ls=":",
            alpha=0.25,
            label="half width",
        )
        axvline(
            x=frq_center + frq_half,
            color="k",
            ls=":",
            alpha=0.25,
            label="half width",
        )
    return frq_center, frq_half


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
    r"""Determine the center of the echo via hermitian symmetry of the time
    domain.

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
    enable_refinement: boolean
        Do not use -- an attempt to correct for cases where we don't get
        the lowest residual out of the Fourier technique.
        The idea was that this was based on the fact that the correlation
        function is calculated between the two upsampled (sinc
        interpolated) functions, and so we should explicitly calculate
        the residual in the original resolution.
        But, there appears to be an error that would need to be debugged
        -- fails some of the basic tests, and the residual does not
        always appear to be the lowest.

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
        s_timedom.getaxis(direct)[0] == 0.0
    ), """In order to
    calculate the signal energy term correctly, the
    signal must start at t=0  so set the start of the
    acquisition in the *non-aliased* time domain to 0 (something like
    data['t2'] -= acqstart) to avoid confusion"""
    t_dw = s_timedom.get_ft_prop(direct, "dt")
    orig_bounds = s_timedom.getaxis(direct)[r_[0, -1]]
    plot_bounds = (
        orig_bounds  # allow this to be set to ± inf if I want to see all
    )
    # }}}
    # {{{ force the axis to *start* at 0
    #     since we're doing FT here, also extend to
    #     twice the length!
    logging.debug(strm("initial size", ndshape(s_ext), "direct is", direct))
    s_ext.ft(
        direct,
        pad=2 ** int(np.ceil(np.log(ndshape(s_ext)[direct]) / np.log(2)) + 1),
    )
    assert s_ext.get_ft_prop(direct, "start_time") == 0, (
        "FT start point should also be equal to zero -- doing otherwise"
        " doesn't make sense"
    )
    # }}}
    # {{{ move into over-sampled time domain -- do this
    # *after* previous to avoid aliasing glitch
    tukeyfilter = s_ext.fromaxis(direct).run(lambda x: sci_win.tukey(len(x)))
    s_ext *= tukeyfilter
    s_ext.ift(
        direct,
        pad=2
        ** int(
            np.ceil(
                np.log(ndshape(s_timedom)[direct] * upsampling) / np.log(2)
            )
        ),
    )
    s_ext[direct : (orig_bounds[-1], None)] = (
        0  # explicitly zero, in case there are aliased negative times!
    )
    # }}}
    # {{{ now I need to throw out the initial, aliased
    #     portion of the signal -- do this manually by
    #     index -- this is more lines than needed in
    #     order to be explanatory.  Note that I
    #     plot the signal before I actually throw stuff
    #     out
    t_dwos = s_ext.get_ft_prop(direct, "dt")  # oversampled dwell
    min_echo = aliasing_slop * t_dw
    logging.debug(
        strm(
            "aliasing slop is",
            aliasing_slop,
            "t_dw is",
            t_dw,
            "(SW of",
            1 / t_dw / 1e3,
            "kHz) min echo is",
            min_echo,
        )
    )
    min_echo_idx = int(min_echo / t_dwos + 0.5)
    min_echo = min_echo_idx * t_dwos
    logging.debug(
        strm(
            "t_dwos is",
            t_dwos,
            "and rounding to the nearest oversampled dwell, min echo is",
            min_echo,
        )
    )
    if fl is not None:
        fl.push_marker()
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
        forplot[direct] -= min_echo  # so zero is the
        #                             first part of the
        #                             echo after the
        #                             aliasing slop, as
        #                             indicated below
        fl.plot(
            forplot,
            label="echo envelope",
        )
        left_bound = forplot[direct][0]
    s_ext[direct, :-min_echo_idx] = s_ext[direct, min_echo_idx:]
    if fl is not None:
        fl.next("power terms")
        forplot = abs(s_ext).mean_all_but(direct)[direct:plot_bounds]
        fl.plot(
            forplot,
            label="echo envelope",
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
    s_correl.run(lambda x: x**2)
    s_correl.ift(direct)
    s_correl.mean_all_but(direct).run(abs)
    s_correl *= normalization_term
    # }}}
    # {{{ calculate the cost function and determine where the center of the
    # echo is!
    cost_func = abs(
        s_energy - s_correl
    )  # b/c this should not be less than 0, so penalize for numerical error
    # when it's not!
    reasonable_energy_range = s_energy.contiguous(
        lambda x: abs(x) > energy_threshold * abs(x.data).max()
    )[0, :]
    _, reasonable_energy_range[1] = s_energy.contiguous(
        lambda x: abs(x)
        > energy_threshold * energy_threshold_lower * abs(x.data).max()
    )[0, :]
    cost_func = cost_func[direct:reasonable_energy_range]
    cost_func.run(lambda x: x / sqrt(abs(x)))  # based on what we'd seen
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
            color="violet",
            alpha=0.5,
        )
        right_bound = forplot[direct][-1]
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
        # {{{ this makes sure we're zoomed on a decent region for the
        #     Hermitian cost function plot
        plt.gca().set_xlim(
            tuple(r_[left_bound, right_bound] / det_devisor(fl))
        )
        # }}}
    # }}}
    echo_idx = int((echo_peak + min_echo) / t_dw + 0.5)
    shift_range = 4
    if enable_refinement and echo_idx - aliasing_slop > shift_range + 1:
        s_foropt = s_timedom
        shifts = nddata(t_dw * r_[-shift_range:shift_range:300j], "echo shift")
        s_foropt.ft(direct)
        # positive pushes time domain to left --
        # should represent that echo occurred later
        # than expected, and so should be added to the echo_idx*t_dw
        # estimate of the echo time
        s_foropt *= np.exp(1j * 2 * np.pi * shifts * s_foropt.fromaxis(direct))
        s_foropt.ift(direct)
        s_foropt = s_foropt[
            direct,
            aliasing_slop : aliasing_slop + 2 * (echo_idx - aliasing_slop) + 1,
        ]
        # the center is now at echo_idx-aliasing_slop
        # {{{ phasing must be done independently for each echo shift
        ph0 = s_foropt[direct, echo_idx - aliasing_slop]
        ph0 /= abs(ph0)
        # }}}
        s_foropt /= ph0
        s_foropt = s_foropt[direct, shift_range:-shift_range]
        s_foropt -= s_foropt.C[direct, ::-1].run(np.conj)
        s_foropt.run(lambda x: abs(x) ** 2)
        s_foropt.mean_all_but("echo shift")
        s_foropt.run(sqrt)
        if fl is not None:
            fl.next("refinement")
            fl.plot(s_foropt, human_units=False)
        s_foropt = s_foropt.argmin("echo shift").item()
        if fl is not None:
            fl.next("refinement")
            axvline(x=s_foropt)
            axvline(
                x=(echo_peak + min_echo) - t_dw * echo_idx, ls=":", alpha=0.25
            )
            fl.next("power terms")
            axvline(
                x=(t_dw * echo_idx + s_foropt - min_echo) / det_devisor(fl),
                ls="-",
                alpha=0.25,
            )
            fl.pop_marker()
        return t_dw * echo_idx + s_foropt
    else:
        if enable_refinement:
            logging.info(
                "warning: can't do hermitian phasing refinement -- not enough"
                " points"
            )
        if fl is not None:
            fl.pop_marker()
        return echo_peak + min_echo


def determine_sign(
    s, signal_freq_range, direct="t2", signal_pathway=None, fl=None
):
    """Determines the sign of the signal based on the difference between the
    signal with the zeroth order phase correction applied to the entire data
    set vs applying the zeroth order phase correction to each individual
    indirect index.  The sign can be used, so that when multiplied by the
    signal the data has a uniform
    sign along the indirect axis.

    Parameters
    ==========
    s : nddata
        Data with a single (dominant) peak, where you want to return the sign
        of the integral over all the data.
    signal_freq_range :  tuple
        Narrow slice range where signal resides.
    direct : str (default "t2")
        Name of the direct dimension.
    signal_pathway : dict (default None)
        If None, the function will go into the properties of the data looking
        for the "coherence_pathway" property. If that doesn't exist the user
        needs to feed a dictionary containing the coherence transfer pathway
        where the signal resides.

    Returns
    =======
    mysign : nddata
        A dataset with all +1 or -1 (giving the sign of the original signal).
        Does *not* include the `direct` dimension.
    """
    if type(direct) is tuple:
        raise ValueError(
            "direct cannot be a tuple! You are probably using an older version"
            " that required you to pass an indirect dimension"
        )
    # To handle both older data where the coherence pathway was not set as a
    # property
    # and handle newer data where it is set I have the following if/else
    signal_pathway = signal_pathway or s.get_prop("coherence_pathway")
    if signal_pathway is not None:
        if list(signal_pathway.values())[0] in s.dimlabels:
            s = select_pathway(s[direct:signal_freq_range], signal_pathway)
    assert s.get_ft_prop(
        direct
    ), "this only works on data that has been FT'd along the direct dimension"
    s /= zeroth_order_ph(s)
    if fl is not None:
        d_forplot = s.C
    s = s.integrate(direct)
    s /= abs(s)  # now s is the phase of the signal along the indirect
    if fl is not None:
        d_forplot /= s  # individually phased
        fl.next("Individually phase corrected")
        fl.image(d_forplot)
    mysign = s.angle.run(lambda x: np.exp(1j * np.round(x / pi) * pi)).run(
        np.sign  # sign essentially rounds to -1 or +1
    )
    if fl is not None:
        fl.next(
            "phase difference between individually rotated and overall rotated"
        )
        fl.plot(s.angle.set_error(None), "o", label="difference")
        fl.plot(mysign.angle.set_error(None), "o", label="rounded")
    return mysign
