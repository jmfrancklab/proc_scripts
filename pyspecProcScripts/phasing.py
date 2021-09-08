"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
from matplotlib.patches import Ellipse
from scipy.optimize import minimize
from pylab import xlim, subplots, axvline, ylim, sca
import numpy as np
from scipy import linalg
import logging


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
            rotation_vector[0],
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
    frq_range=(-2.5e3, 2.5e3),  # assumes signal is somewhat on resonance
    aliasing_slop=3,
    band_mask=False,
    fl=None,
):
    r"""determine the center of the echo via hermitian symmetry of the time domain.

    Parameters
    ==========
    direct:             str
        Axis of data (i.e., direct dimension).
    frq_range:          tuple of floats
        range over which to slice data in the
        frequency domain for frequency filtering
    aliasing_slop:          int
        Because we sinc interpolate here, we need to allow for the fact that
        the very beginning and ending of the time-domain signal are
        interpolated to match.
        In units of the dwell time of the original signal, this slices out the
        aliased parts at the beginning and end of the time-domain signal.
    band_mask:          boolean
        determines the type of mask used on the 2D
        residual for accurate calculation of cost
        function. The mask is used to circumvent
        aliasing of the cost function.
        When True, applies rectangular mask to the
        2D residual.
        When False, applies triangular mask to the
        2D residual.
    band_mask_no:       int
        determines number of points for rectangular mask
    final_range:        tuple
        range that will slice out everything but the very first few points
        which gives an artificial minimum and messes with the
        cost function.
    """
    orig_dt = s.get_ft_prop(direct, "dt")
    s.ft(direct)
    s = s[direct:frq_range]
    s.ift(direct, pad=1024 * 16)
    new_dt = s.get_ft_prop(direct, "dt")
    non_aliased_range = r_[aliasing_slop,-aliasing_slop]*int(orig_dt/new_dt)
    ini_start = s.getaxis(direct)[0]
    s = s[direct,non_aliased_range[0]:non_aliased_range[1]]
    ini_delay = s.getaxis(direct)[0]
    logger.debug(strm("ini delay is",ini_delay))
    if fl is not None:
        fl.push_marker()
        fig, ax_list = subplots(2, 3, figsize=(15, 15))
        fl.next("Hermitian Function Test Diagnostics", fig=fig)
        fl.plot(abs(s), ax=ax_list[0, 0], human_units=False)
        ax_list[0, 0].set_title("Data with Padding")
    probable_center = abs(s).convolve(direct,orig_dt*3).argmax(direct).item() # convolve just for some signal averaging
    residual = s[direct:(ini_start,ini_start+(probable_center-ini_start)*2)]
    if fl is not None:
        fl.plot(abs(residual), ':', ax=ax_list[0, 0], human_units=False)
    N = ndshape(residual)[direct]
    mid_idx = N // 2 + N % 2 - 1
    residual = residual[direct, 0 : 2 * mid_idx + 1]
    if fl is not None:
        if band_mask:
            title_str = "rectangular mask"
        else:
            title_str = "triangular mask"
        fl.next("cost function %s - freq filter" % title_str)
        residual.name("absolute value")
        fl.plot(
            abs(residual)
            .mean_all_but(direct)
            .rename(direct, "center")
            .set_units("center", "s"),
            color="k",
            alpha=0.5,
            human_units=False,
        )
    dt = residual.get_ft_prop(direct, "dt")
    # the shifts themselves run from 0 to mid_idx -- the echo-centers
    # these correspond to are different. Also, we're not really
    # worried about the sub-integer shift here, because we can just
    # use sinc interpolation before we start -- see:
    # https://jmfrancklab.slack.com/archives/CLMMYDD98/p1623354066039100
    shifts = nddata(dt * (r_[0:mid_idx]), "shift")
    shifts.set_units("shift", "s")
    logger.debug(strm("Length of shifts dimension:", ndshape(shifts)["shift"]))
    residual.ft(direct)
    residual *= np.exp(-1j * 2 * pi * shifts * residual.fromaxis(direct))
    residual.ift(direct)
    logger.debug(strm("Length of t2 dimension:", ndshape(residual)[direct]))
    assert (
        ndshape(residual)[direct] % 2 == 1
    ), "t2 dimension *must* be odd, please check what went wrong."
    # {{{phase correct 
    center_point = residual[direct, mid_idx]
    residual /= center_point/abs(center_point)
    # }}}
    if fl is not None:
        if "power" in s.dimlabels:
            fl.image(residual["power", -4], ax=ax_list[0, 1])
        else:
            fl.image(residual, ax=ax_list[0, 1])
        ax_list[0, 1].set_title("shifted and phased")
    # I tried removing the following line, and the 2D residual actually
    # looks sharper, but then the 1D plot of the cost function disappears
    # -- what's up with that?
    residual = abs(residual - residual[direct, ::-1].runcopy(np.conj)).mean_all_but(
        ["shift", direct]
    )
    if fl is not None:
        fl.image(residual, ax=ax_list[0, 2])
        ax_list[0, 2].set_title("Residual 2D")
    if band_mask:
        n_mask_pts = int(band_mask / dt)
        if n_mask_pts % 2:
            n_mask_pts -= 1
        residual[direct, mid_idx + n_mask_pts :] = 0
        residual[direct, : mid_idx - n_mask_pts] = 0
    else:
        # Here we would do the calculation outlined for the triangle mask,
        # but with only two new variables -- A and B
        A = r_[-mid_idx : mid_idx + 1]
        A = -abs(A)
        A += mid_idx
        B = residual.fromaxis("shift").data / dt
        A = A.reshape(-1, 1)
        B = B.reshape(1, -1)
        mask = np.greater_equal(A, B)
        mask = np.multiply(mask, 1)
        mask = mask.astype(float)
        norm = np.floor(mid_idx - B)
        norm = norm.astype(float)
        mask /= ((norm) * 2)**2
        mask = nddata(mask, [direct, "shift"])
        mask.setaxis(direct, residual.getaxis(direct))
        mask.setaxis("shift", residual.getaxis("shift"))
        if fl is not None:
            fl.image(mask, ax=ax_list[1, 0], human_units=False)
            ax_list[1, 0].set_title("Mask")
        residual *= mask
    if fl is not None:
        fl.image(residual, ax=ax_list[1, 1])
        ax_list[1, 1].set_title("Masked Residual")
    residual.rename("shift", "center").setaxis(
        "center", dt * (mid_idx - r_[0:mid_idx]) + ini_delay
    )
    residual = residual["center", ::-1]
    if fl is not None:
        fl.image(residual, ax=ax_list[1, 2], human_units=False)
        ax_list[1, 2].set_title("Masked Residual -- Relabeled")
        fig.tight_layout()
    residual.mean(direct)
    residual.set_units("center", "s")
    best_shift = residual.C.argmin("center").item()
    # slices first few points out as usually the artifacts give a minimum
    if fl is not None:
        fl.next("cost function %s - freq filter" % title_str)
        fl.twinx(orig=False, color="red")
        residual.name("cost function")
        fl.plot(residual, color="r", alpha=0.5, human_units=False)
        fl.plot(residual["center" : (best_shift - 4e-3, best_shift)], ':', alpha=0.5, human_units=False)
        #ylim(0, residual["center" : (best_shift - 4e-3, best_shift)].data.max())
        print("I find max",residual["center" : (best_shift - 4e-3, best_shift)].data.max())
        axvline(x=best_shift, c="k", linestyle="--")
        for dwell_int in r_[-5:5]:
            axvline(
                x=best_shift - (orig_dt * dwell_int), alpha=0.4, c="k", linestyle=":"
            )
        fl.pop_marker()
    return best_shift
