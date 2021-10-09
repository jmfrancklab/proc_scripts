"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
from matplotlib.patches import Ellipse
from matplotlib.transforms import blended_transform_factory
from scipy.optimize import minimize
from pylab import xlim, subplots, axvline, ylim, sca
import numpy as np
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
    aliasing_slop=3,
    band_mask=False,
    fl=None,
):
    r"""determine the center of the echo via hermitian symmetry of the time domain.

    Parameters
    ==========
    direct:             str
        Axis of data (i.e., direct dimension).
    aliasing_slop:          int
        Because we sinc interpolate here, we need to allow for the fact that
        the very beginning and ending of the time-domain signal are
        interpolated to match.
        This value is the multiple of the dwell time of the original signal
        that is sliced out due to the fact this amount of signal is aliased 
        at the beginning and end of the time-domain signal.
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
    s = s.C # work on a copy, so that we're not messing w/ anything
    orig_dt = s.get_ft_prop(direct, "dt")
    if not s.get_ft_prop(direct):
        s.ft(direct)
    s.ift(direct, pad=ndshape(s)[direct]*50)
    new_dt = s.get_ft_prop(direct, "dt")
    non_aliased_range = r_[aliasing_slop,-aliasing_slop]*int(orig_dt/new_dt)
    if aliasing_slop > 0:
        s = s[direct,non_aliased_range[0]:non_aliased_range[1]]
    s_envelope = s.C.mean_all_but(direct).run(abs)
    if fl is not None:
        fl.push_marker()
        fig, ax_list = subplots(2, 2, figsize=(15, 15))
        fl.next("Hermitian Function Test Diagnostics", fig=fig)
        fl.plot(s_envelope, ax=ax_list[0, 0], human_units=False)
        ax_list[0, 0].set_title("Data with Padding")
    #s_envelope.convolve(direct, orig_dt*3)
    #if fl is not None:
    #    fl.plot(s_envelope, ax=ax_list[0, 0], human_units=False)
    #peak_triple = abs(s).mean_all_but(direct).convolve(direct,orig_dt*3).contiguous(lambda x: x > 0.25 * x.data.max())[0,:]
    peak_triple = s_envelope.contiguous(lambda x: x > 0.25 * x.data.max())[0,:]
    # {{{ the peak triple gives the left and right thresholds, with the
    #     peak max in the middle
    #     we move the left point as needed to make sure that peak is in
    #     the first half of the resulting slice
    peak_triple = r_[ peak_triple[0],
            s[direct:peak_triple].mean_all_but(direct).run(abs).argmax().item(),
            peak_triple[1] ]
    if fl is not None:
        for j in range(3):
            ax_list[0, 0].axvline(x=peak_triple[j], ls=':', color='k', alpha=0.5, linewidth=2)
    extension_left, extension_right = abs(np.diff(peak_triple))
    if extension_right < extension_left + orig_dt:
        if peak_triple[2] < s.getaxis(direct)[-1]:
            peak_triple[2] = s.getaxis(direct)[-1]
            extension_left, extension_right = abs(np.diff(peak_triple))
        if extension_right < extension_left + orig_dt:
            extension_left = extension_right - orig_dt
            peak_triple[0] = peak_triple[1] - extension_left
    s = s[direct:peak_triple[r_[0,-1]]]
    # }}}
    if fl is not None:
        ax_list[0,0].axvspan(0,orig_dt,color='k',alpha=0.1)
        for j in range(3):
            ax_list[0, 0].axvline(x=peak_triple[j], color='k', alpha=0.5, linewidth=2)
        fl.plot(abs(s), ':', ax=ax_list[0, 0], linewidth=4, human_units=False)
    direct_startpoint = s.getaxis(direct)[0]
    logger.debug(strm("ini delay is",direct_startpoint))
    N = ndshape(s)[direct]
    mid_idx = N // 2 + N % 2 - 1
    s = s[direct, 0 : 2 * mid_idx + 1]
    if fl is not None:
        if band_mask:
            title_str = "rectangular mask"
        else:
            title_str = "triangular mask"
        fl.next("cost function %s - freq filter" % title_str)
        s.name("absolute value")
        fl.plot(
            abs(s)
            .mean_all_but(direct)
            .rename(direct, "center")
            .set_units("center", "s"),
            color="k",
            alpha=0.5,
            human_units=False,
        )
    dt = s.get_ft_prop(direct, "dt")
    # the shifts themselves run from 0 to mid_idx -- the echo-centers
    # these correspond to are different. Also, we're not really
    # worried about the sub-integer shift here, because we can just
    # use sinc interpolation before we start -- see:
    # https://jmfrancklab.slack.com/archives/CLMMYDD98/p1623354066039100
    shifts = nddata(dt * (r_[0:mid_idx]), "shift")
    shifts.set_units("shift", "s")
    logger.debug(strm("Length of shifts dimension:", ndshape(shifts)["shift"]))
    s.ft(direct)
    s *= np.exp(-1j * 2 * pi * shifts * s.fromaxis(direct))
    s.ift(direct)
    logger.debug(strm("Length of t2 dimension:", ndshape(s)[direct]))
    assert (
        ndshape(s)[direct] % 2 == 1
    ), "t2 dimension *must* be odd, please check what went wrong."
    # {{{phase correct 
    center_point = s[direct, mid_idx]
    s /= center_point/abs(center_point)
    # }}}
    if fl is not None:
        if "power" in s.dimlabels:
            fl.image(s["power", -4], ax=ax_list[0, 1])
        else:
            fl.image(s, ax=ax_list[0, 1])
        ax_list[0, 1].set_title("shifted and phased")
    s = abs(s - s[direct, ::-1].runcopy(np.conj)).mean_all_but(
        ["shift", direct]
    )
    if fl is not None:
        fl.image(s, ax=ax_list[0, 1])
        ax_list[0, 1].set_title("Residual 2D")
    # {{{ the "center axis makes more sense, so implement it sooner"
    shift_idxunits = s.fromaxis("shift").data / dt # used below
    #s /= abs(center_point.C.mean_all_but(["shift",direct])) # weight the cost by the amplitude of the center point
    s.rename("shift", "center").setaxis(
        "center", dt * (mid_idx - r_[0:mid_idx]) + direct_startpoint
    )
    s.set_units("center", "s")
    if fl is not None:
        fl.image(s.C.cropped_log()['center',::-1], ax=ax_list[1, 0])
        ax_list[1, 0].set_title("Residual 2D\nrelabel, signal weight, and cropped log")
    # }}}
    if band_mask:
        n_mask_pts = int(band_mask / dt)
        if n_mask_pts % 2:
            n_mask_pts -= 1
        s[direct, mid_idx + n_mask_pts :] = 0
        s[direct, : mid_idx - n_mask_pts] = 0
    else:
        # Here we would do the calculation outlined for the triangle mask,
        # but with only two new variables -- A and B
        A = r_[-mid_idx : mid_idx + 1]
        A = -abs(A)
        A += mid_idx
        B = shift_idxunits 
        A = A.reshape(-1, 1)
        B = B.reshape(1, -1)
        mask = np.greater_equal(A, B)
        mask = np.multiply(mask, 1)
        mask = mask.astype(float)
        norm = np.floor(mid_idx - B)
        norm = norm.astype(float)
        mask /= ((norm) * 2)**2
        mask = nddata(mask, [direct, "center"])
        mask.setaxis(direct, s.getaxis(direct))
        mask.setaxis("center", s.getaxis("center"))
        s *= mask
    if fl is not None:
        fl.image(s.C.cropped_log()["center",::-1], ax=ax_list[1, 1])
        ax_list[1, 1].set_title("cropped log of Masked Residual")
        fig.tight_layout()
    s = s["center", ::-1]
    s.mean_all_but(["center"])
    best_shift = s.C.argmin("center").item()
    # slices first few points out as usually the artifacts give a minimum
    if fl is not None:
        fl.next("cost function %s - freq filter" % title_str)
        fl.twinx(orig=False, color="red")
        s.name("cost function")
        fl.plot(s, color="r", alpha=0.5, human_units=False)
        cost_scale = s['center':(-4*orig_dt+best_shift,4*orig_dt+best_shift)].max().item()
        ylim((0,cost_scale))
        axvline(x=best_shift, c="r", linestyle="--")
        ax = plt.gca()
        trans = blended_transform_factory(x_transform=ax.transData, y_transform=ax.transAxes)
        text(x = best_shift+0.0005,
                y = 0.8,
                s = "best shift is: %f"%best_shift,
                fontsize = 16,
                color='red',
                transform=trans,
                )
        for dwell_int in r_[-5:5]:
            axvline(
                x=best_shift - (orig_dt * dwell_int), alpha=0.4, c="r", linestyle=":"
            )
        fl.pop_marker()
    return best_shift
