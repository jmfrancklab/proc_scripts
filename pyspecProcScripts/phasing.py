"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
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
    aliasing_slop=3,
    band_mask=False,
    amp_threshold=0.5,
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
        fig_forlist, ax_list = plt.subplots(2, 2, figsize=(15, 15))
        fl.next("Hermitian Function Test Diagnostics")
        fig_forlist.suptitle(" ".join(["Hermitian Diagnostic"] + [j for j in [fl.basename] if j is not None]))
        fl.plot(s_envelope, ax=ax_list[0, 0], human_units=False)
        ax_list[0, 0].set_title("Data with Padding")
    peak_triple = s_envelope.contiguous(lambda x: x > amp_threshold * x.data.max())[0,:]
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
    # {{{ make sure we don't test for center too close to either edge
    if peak_triple[0] < s.getaxis(direct)[0]+orig_dt:
        peak_triple[0] = s.getaxis(direct)[0]+orig_dt
    if peak_triple[2] > s.getaxis(direct)[-1]-orig_dt:
        peak_triple[2] = s.getaxis(direct)[-1]-orig_dt
    # }}}
    slice_start,slice_stop = s.get_range(direct,*peak_triple[r_[0,-1]])
    # }}}
    if fl is not None:
        ax_list[0,0].axvspan(0,orig_dt,color='k',alpha=0.1)
        for j in range(3):
            ax_list[0, 0].axvline(x=peak_triple[j], color='k', alpha=0.5, linewidth=2)
        fl.plot(abs(s)[direct,slice(slice_start,slice_stop)], ':', ax=ax_list[0, 0], linewidth=4, human_units=False)
    N = ndshape(s)[direct]
    direct_startpoint = s.getaxis(direct)[0]
    center_startpoint = s.getaxis(direct)[slice_start]
    if fl is not None:
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
    # {{{ put test_center at the start/bottom of the data
    test_center = nddata(s.getaxis(direct)[slice_start:slice_stop],'center')
    s.rename(direct,'offset')
    s.setaxis('offset', lambda x: x-direct_startpoint) # so it starts at 0 now
    s.ft('offset')
    # the following confused me a little -- I wanted to take the difference between test_center and the start of the axis,
    # but pyspecdata handles issues about the start of the axis,
    # so we just shift so that test_center moves to zero
    s *= np.exp(1j * 2 * pi * test_center * s.fromaxis('offset'))
    s.ift('offset')
    phcorr = s['offset', 0]
    phcorr /= abs(phcorr)
    s /= phcorr
    # }}}
    if fl is not None:
        fl.image(s, ax=ax_list[0, 1])
        ax_list[0, 1].set_title("shifted and phased")
    # {{{ calculate the mask
    center_idx = r_[slice_start:slice_stop] # test these points to see if they are the center
    # Now, excluding the center point, determine how many points our "mask"
    rms_size = N - center_idx - 1 # number of points between (not including) the center and the end of the data
    rms_size[rms_size > center_idx] = center_idx[rms_size > center_idx] # if there are fewer points on the other side, we are limited by those
    mymask = nddata(r_[1:N],'offset') <= nddata(rms_size,'center')
    # }}}
    if band_mask:
        raise ValueError("band mask no longer supported -- use an earlier version")
    s = abs(s['offset',1:] - s['offset', :0:-1].runcopy(np.conj)).mean_all_but(
        ["center", 'offset']
    )
    if fl is not None:
        fl.image(s, ax=ax_list[1, 0])
        ax_list[1, 0].set_title("Residual 2D")
    s *= mymask
    s /= nddata(rms_size,'center')**2 # it's an L1 difference, but this is
    #                                   squared - this is just what worked
    #                                   well -- we could use a rationale for
    #                                   the squared, actually
    #                                   → maybe once for normalization, once for probability?
    if fl is not None:
        fl.image(s.C.cropped_log(), ax=ax_list[1, 1])
        ax_list[1, 1].set_title("cropped log of Masked Residual")
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
    assert s.get_ft_prop(direct), "this only works on data that has been FT'd along the direct dimension"
    if fl is not None:
        fl.push_marker()
        fl.next('selected pathway')
        fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
    data_sgn = s.C.sum(direct)
    data_sgn /= zeroth_order_ph(data_sgn)
    data_sgn.run(np.real).run(lambda x: np.sign(x))
    if fl is not None:
        fl.next('check sign')
        fl.image((s.C.setaxis(
'vd','#').set_units('vd','scan #'))*data_sgn)
        fl.pop_marker()
    return data_sgn
