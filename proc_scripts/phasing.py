"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
from matplotlib.patches import Ellipse
from scipy.optimize import minimize
def zeroth_order_ph(d, fl=None):
    r'''determine the covariance of the datapoints
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
    '''
    cov_mat = cov(c_[
        d.data.real.ravel(),
        d.data.imag.ravel()].T,
        aweights=abs(d.data).ravel()**2 # when running proc_square_refl, having
        # this line dramatically reduces the size of the imaginary component
        # during the "blips," and the magnitude squared seems to perform
        # slightly better than the magnitude -- this should be both a robust
        # test and a resolution for issue #23
        )
    eigenValues, eigenVectors = eig(cov_mat)
    mean_point = d.data.ravel().mean()
    mean_vec = r_[mean_point.real,mean_point.imag]
    # next 3 lines from stackexchange -- sort by
    # eigenvalue
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx] # first dimension x,y second evec #
    # determine the phase angle from direction of the
    # largest principle axis plus the mean
    # the vector in the direction of the largest
    # principle axis would have a norm equal to the
    # sqrt (std not variance) of the eigenvalue, except
    # that we only want to rotate when the distribution
    # is assymetric, so include only the excess of the
    # larger eval over the smaller
    assymetry_mag = sqrt(eigenValues[0])-sqrt(eigenValues[1])
    try:
        assym_ineq = (assymetry_mag*eigenVectors[:,0]*mean_vec).sum()
    except:
        raise ValueError(strm("check the sizes of the following:",size(assymetry_mag),size(eigenVectors),size(mean_vec)))
    if assym_ineq > 0:
        # we want the eigenvector on the far side of the ellipse
        rotation_vector = mean_vec + assymetry_mag*eigenVectors[:,0]
    else:
        rotation_vector = mean_vec - assymetry_mag*eigenVectors[:,0]
    ph0 = arctan2(rotation_vector[1],rotation_vector[0])
    if fl:
        d_forplot = d.C
        fl.next('check covariance test')
        fl.plot(
                d_forplot.data.ravel().real,
                d_forplot.data.ravel().imag,
                '.',
                alpha=0.25,
                label='before'
                )
        d_forplot /= exp(1j*ph0)
        fl.plot(
                d_forplot.data.ravel().real,
                d_forplot.data.ravel().imag,
                '.',
                alpha=0.25,
                label='after'
                )
        fl.plot(0,0,'ko', alpha=0.5)
        fl.plot(mean_vec[0],mean_vec[1],'kx', label='mean', alpha=0.5)
        evec_forplot = sqrt(eigenValues.reshape(1,2))*ones((2,1))*eigenVectors # scale by the std, not the variance!
        evec_forplot += mean_vec.reshape((-1,1))*ones((1,2))
        fl.plot(evec_forplot[0,0],evec_forplot[1,0],'o', alpha=0.5,
                label='first evec')
        fl.plot(evec_forplot[0,1],evec_forplot[1,1],'o', alpha=0.5)
        fl.plot(rotation_vector[0],rotation_vector[1],'o', alpha=0.5,
                label='rotation vector')
        norms = sqrt((evec_forplot**2).sum(axis=0))
        ell = Ellipse(xy=mean_vec,
                width=2*sqrt(eigenValues[0]),
                height=2*sqrt(eigenValues[1]),
                angle=180/pi*arctan2(eigenVectors[1,0],
                    eigenVectors[0,0]),
                color='k', fill=False)
        ax = gca()
        ax.set_aspect('equal', adjustable='box')
        ax.add_patch(ell)
    return exp(1j*ph0)

def ph1_real_Abs(s,dw):
    r''' Performs first order phase correction with cost function
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
    '''
    fl.push_marker() 
    ph1 = nddata(r_[-5:5:70j]*dw,'phcorr')
    dx = diff(ph1.getaxis('phcorr')[r_[0,1]]).item()
    ph1 = exp(-1j*2*pi*ph1*s.fromaxis('t2'))
    s_cost = s * ph1
    ph0 = s_cost.C.sum('t2')
    ph0 /= abs(ph0)
    s_cost /= ph0
    fl.next('phasing cost function')
    s_cost.run(real).run(abs).sum('t2')
    fl.plot(s_cost,'.')
    ph1_opt = s_cost.argmin('phcorr').item()
    print('optimal phase correction',repr(ph1_opt))
    # }}}
    # {{{ apply the phase corrections
    def applyphase(arg,ph1):
        arg *= exp(-1j*2*pi*ph1*arg.fromaxis('t2'))
        ph0 = arg.C.sum('t2')
        ph0 /= abs(ph0)
        arg /= ph0
        return arg
    def costfun(ph1):
        if type(ph1) is ndarray:
            ph1 = ph1.item()
        temp = s.C
        retval = applyphase(temp,ph1).run(real).run(abs).sum('t2').item()
        return retval
    print("rough opt cost function is",costfun(ph1_opt))
    r = minimize(costfun,ph1_opt,
            bounds=[(ph1_opt-dx,ph1_opt+dx)])
    assert r.success
    s = applyphase(s,r.x.item())
    fl.plot(r.x,r.fun,'x')
    fl.pop_marker()
    return s

def hermitian_function_test(s, down_from_max=0.5, shift_val=1.0, fl=None):

    r"""Determine the center of the echo by assuming that the echo has Hermitian symmetry.

    As part of the algorithm, the center point must be set to entirely real.  This is done using :func:`zeroth_order_ph`

    .. todo::
        Alex: please update this using standard Parameters and
        Returns blocks
    """
    if s.get_ft_prop('t2', ['start','freq']) is None:
        # this sets the frequency startpoint appropriately for an FT --> it
        # should be possible to do this w/out actually performing an ft
        s.ft('t2', shift=True)
        s.ift('t2')
    # {{{ determine where the "peak" of the echo is,
    # and use it to determine the max
    # shift
    s = s.C # need to copy, since I'm manipulating the axis here,
    # and am assuming the axis of the source data is not manipulated
    data_for_peak = abs(s).mean_all_but(['t2'])
    max_val = data_for_peak.data.max()
    pairs = data_for_peak.contiguous(lambda x: abs(x) >
            max_val*down_from_max)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    peak_center = peak_location.mean()
    s.setaxis('t2',lambda x: x-peak_center)
    s.register_axis({'t2':0})
    max_shift = diff(peak_location).item()/2
    # }}}
    # {{{ construct test arrays for T2 decay and shift
    shift_t = nddata(r_[-1*shift_val:1*shift_val:1200j]*max_shift, 'shift')
    # }}}
    s_foropt = s.C
    # {{{ time shift and correct for T2 decay
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*
            s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    # }}}
    # {{{ make sure there's and odd number of points
    # and set phase of center point to 0
    logger.info(strm("here is t2",s_foropt.getaxis('t2')))
    logger.info(strm(ndshape(s_foropt)))
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
    logger.info(strm(ndshape(s_foropt)))
    n_points = ndshape(s_foropt)['t2']
    if n_points % 2 == 0:
        s_foropt = s_foropt['t2',:-1]
        n_points -= 1
    logger.info(strm(ndshape(s_foropt)))
    center_point = s_foropt['t2',n_points//2+1]
    s_foropt /= zeroth_order_ph(center_point)
    # }}}
    residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).mean_all_but(['shift','R2'])
    # in the following, weight for the total signal recovered
    residual = residual / abs(center_point).mean_all_but(['shift','R2'])
    residual.reorder('shift')
    minpoint = residual.argmin()
    best_shift = minpoint['shift']

    if fl:
        fl.next('residual for hermitian test')
        fl.plot(residual)
    return best_shift+peak_center

