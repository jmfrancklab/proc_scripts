"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
from matplotlib.patches import Ellipse
from scipy.optimize import minimize
from pylab import xlim
import numpy as np
import pyspecdata as pysp
from scipy import linalg
import logging
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
    cov_mat = np.cov(c_[
        d.data.real.ravel(),
        d.data.imag.ravel()].T,
        aweights=abs(d.data).ravel()**2 # when running proc_square_refl, having
        # this line dramatically reduces the size of the imaginary component
        # during the "blips," and the magnitude squared seems to perform
        # slightly better than the magnitude -- this should be both a robust
        # test and a resolution for issue #23
        )
    eigenValues, eigenVectors = linalg.eig(cov_mat)
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
    assymetry_mag = float(sqrt(eigenValues[0])-sqrt(eigenValues[1]))
    try:
        assym_ineq = (assymetry_mag*eigenVectors[:,0]*mean_vec).sum()
    except:
        raise ValueError(strm("check the sizes of the following:",size(assymetry_mag),size(eigenVectors),size(mean_vec)))
    if assym_ineq > 0:
        # we want the eigenvector on the far side of the ellipse
        rotation_vector = mean_vec + assymetry_mag*eigenVectors[:,0]
    else:
        rotation_vector = mean_vec - assymetry_mag*eigenVectors[:,0]
    ph0 = np.arctan2(rotation_vector[1],rotation_vector[0])
    if fl is not None:
        d_forplot = d.C
        fl.next('check covariance test')
        fl.plot(
                d_forplot.data.ravel().real,
                d_forplot.data.ravel().imag,
                '.',
                alpha=0.25,
                label='before'
                )
        d_forplot /= np.exp(1j*ph0)
        fl.plot(
                d_forplot.data.ravel().real,
                d_forplot.data.ravel().imag,
                '.',
                alpha=0.25,
                label='after'
                )
        fl.plot(0,0,'ko', alpha=0.5)
        fl.plot(mean_vec[0],mean_vec[1],'kx', label='mean', alpha=0.5)
        evec_forplot = sqrt(eigenValues.reshape(1,2))*np.ones((2,1))*eigenVectors # scale by the std, not the variance!
        evec_forplot += mean_vec.reshape((-1,1))*np.ones((1,2))
        fl.plot(evec_forplot[0,0],evec_forplot[1,0],'o', alpha=0.5,
                label='first evec')
        fl.plot(evec_forplot[0,1],evec_forplot[1,1],'o', alpha=0.5)
        fl.plot(rotation_vector[0],rotation_vector[1],'o', alpha=0.5,
                label='rotation vector')
        norms = sqrt((evec_forplot**2).sum(axis=0))
        ell = Ellipse(xy=mean_vec,
                width=2*sqrt(eigenValues[0]),
                height=2*sqrt(eigenValues[1]),
                angle=180/pi*np.arctan2(eigenVectors[1,0],
                    eigenVectors[0,0]),
                color='k', fill=False)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        ax.add_patch(ell)
    return np.exp(1j*ph0)
def ph1_real_Abs(s,dw,ph1_sel=0,ph2_sel=1,fl = None):
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
    if fl: fl.push_marker() 
    ph1 = nddata(r_[-5:5:70j]*dw,'phcorr')
    dx = np.diff(ph1.getaxis('phcorr')[r_[0,1]]).item()
    ph1 = np.exp(-1j*2*pi*ph1*s.fromaxis('t2'))
    s_cost = s * ph1
    ph0 = s_cost.C.sum('t2')
    ph0 /= abs(ph0)
    s_cost /= ph0
    if fl:
        fl.next('phasing cost function')
    s_cost.run(np.real).run(abs).sum('t2')
    s_cost = s_cost['ph2',1]['ph1',0]
    if fl:
        fl.plot(s_cost,'.')
    s_cost.smoosh(['vd','phcorr'],'transients')
    ph1_opt = np.asarray(s_cost.argmin('transients').item())
    logging.info(strm("THIS IS PH1_OPT",ph1_opt))
    logging.info(strm('optimal phase correction',repr(ph1_opt)))
    # }}}
    # {{{ apply the phase corrections
    def applyphase(arg,ph1):
        arg *= np.exp(-1j*2*pi*ph1*arg.fromaxis('t2'))
        ph0 = arg.C.sum('t2')
        ph0 /= abs(ph0)
        arg /= ph0
        return arg
    def costfun(ph1):
        logging.info(strm("PH1 IS THIS:",ph1))
        if type(ph1) is np.ndarray:
            ph1 = ph1.item()
        temp = s.C
        retval = applyphase(temp,ph1).run(np.real).run(abs).sum('t2').item()
        return retval
    r = minimize(costfun,
            ph1_opt,
            bounds = [(ph1_opt-dx,ph1_opt+dx)])
    assert r.success
    s = applyphase(s,r.x.item())
    fl.plot(r.x,r.fun,'x')
    fl.pop_marker()
    return s
    #}}}
def hermitian_function_test(s, down_from_max=0.5, rel_shift=1.5,shift_points=1200, fl=None):
    r"""determine the center of the echo

    Parameters
    ==========
    down_from_max:  float
        Slice determine the portion of the time axis where the echo
        is at this fraction of the echo max, and use this as the 
        "peak bounds" of the echo.

        Use the width of this peak location as the width of the window
        used for the Hermitian function test calculation
   rel_shift:       float
        When testing the cost function, shift up to half the time
        range (see down_from_max) in either direction, multiplied
        (expanded/contracted) by this number
    """
    #{{{ determine where the "peak" of the echo is,
    # and use it to determine the max shift
    s = s.C # need to copy, since I'm manipulating the axis here
    dt = np.diff(s.getaxis('t2')[:2]).item()
    logging.debug(strm("For hermitian test, my dwell time is",dt))
    # am assuming the axis of the source data is not manipulated
    data_for_peak = abs(s).mean_all_but(['t2'])
    max_val = data_for_peak.data.max()
    data_for_peak /= max_val
    pairs = data_for_peak.contiguous(lambda x: abs(x) >
            down_from_max)
    peak_bounds = pairs[0,:] #gives the longest pair
    peak_center = data_for_peak['t2':(0.1e-3,None)].argmax('t2').item()
    s.setaxis('t2',lambda x: x-peak_center)
    s.register_axis({'t2':0})
    logging.debug(strm("closest to 0", s.getaxis('t2')[np.argmin(abs(s.getaxis('t2')-0))]))
    max_shift = np.min(abs(peak_center-peak_bounds))
    #}}}
    #{{{construct test arrays for T2 decay and shift
    if np.isscalar(rel_shift):
        rel_shift = [-rel_shift,rel_shift]
    shift_t = nddata(r_[rel_shift[0]:rel_shift[1]:shift_points*1j]*max_shift,'shift')
    data_for_peak.setaxis('t2',lambda x:x-peak_center) #for plotting later
    logging.debug(strm("closest to 0", s.getaxis('t2')[np.argmin(abs(s.getaxis('t2')-0))]))
    #data_for_peak = data_for_peak['t2':(rel_shift[0]*max_shift,rel_shift[1]*max_shift)]
    #}}}
    #{{{time shift and correct for T2 decay
    orig_t2 = ndshape(s)['t2']
    s.ft('t2',pad=orig_t2*2)
    s *= np.exp(1j*2*pi*shift_t*
            s.fromaxis('t2'))
    s.set_units('shift','s')
    s.ift('t2')
    s = s['t2',:orig_t2]
    logging.debug(strm("closest to 0", s.getaxis('t2')[np.argmin(abs(s.getaxis('t2')-0))]))
    if fl:
        fl.next('shifted data')
        fl.image(s.C.mean_all_but(['shift','t2']))
    #}}}
    #{{{make sure there's an odd number of points and set phase of center point to 0
    logging.debug(strm("before taking the slice, endpoints are",s.getaxis('t2')[r_[0,-1]]))
    logging.debug(strm("closest to 0", s.getaxis('t2')[np.argmin(abs(s.getaxis('t2')-0))]))
    # {{{ this passes the assertion, but I
    #     need jf_hack for the next step to
    #     work -- this means the axis isn't
    #     made of numerically exactly equal
    #     steps, for some reason
    jf_hack = True
    test = np.diff(s.getaxis('t2'))
    assert np.allclose(dt*np.ones_like(test),test)
    if jf_hack:
        x = s.getaxis('t2')
        x[:] /= dt
        x[:] = np.round(x)
        x *= dt
    # }}}
    s_hermitian = s['t2':(-max_shift,max_shift)].C
    logging.debug(strm("closest to 0", s.getaxis('t2')[np.argmin(abs(s.getaxis('t2')-0))]))
    n_points = ndshape(s_hermitian)['t2']
    logging.debug(strm('size of axis',len(s_hermitian.getaxis('t2'))))
    logging.debug(strm('steps',np.unique(np.diff(s_hermitian.getaxis('t2')))))
    logging.debug(strm('full list of axis after slice',s_hermitian.getaxis('t2')))
    logging.debug(strm("check for endpoints",s_hermitian.getaxis('t2')[r_[0,-1]],
        s_hermitian.getaxis('t2')[0] + s_hermitian.getaxis('t2')[-1],
        s_hermitian.getaxis('t2')[n_points//2]))
    if n_points % 2 == 0:
        raise ValueError("Getting an even slice -- shouldn't happen if your "
                "axis coord are in register with t=0 and you take a slice "
                "(-max_shift,max_shift)")
    logging.debug("check for center", n_points, s_hermitian.getaxis('t2')[n_points//2:n_points//2+3])
    assert s_hermitian.getaxis('t2')[n_points//2] == 0
    center_point = s_hermitian['t2',n_points//2]
    logging.info(strm(center_point/abs(center_point)))
    #I want to do this center_point = s_hermitian['t2':0], but why not be consistent
    s_hermitian /= center_point/abs(center_point)
    if fl:
        fl.next('hermitian data')
        fl.image(s_hermitian.C.mean_all_but(['shift','t2']))
    #}}}
    #though we want to average the residual for each FID, try making an average
    #FID and calculating the residual, as I do in the test towards the end
    #mean_before = s_hermitian.C.mean_all_but(['shift','t2'])
    #residual = abs(mean_before - mean_before['t2',::-1].runcopy(conj))
    # 
    # I have something like ΣⱼΣᵢ |aᵢⱼ-bᵢⱼ|² -- (say i vector subscript and j
    # repeats) -- pushing the Σⱼ inside the |...|² does matter, so I do restore
    # the original code here, though it's likely slower for large datasets
    #
    # actually, we find empirically that ΣⱼΣᵢ |aᵢⱼ-bᵢⱼ|,
    # rather than ΣⱼΣᵢ |aᵢⱼ-bᵢⱼ|²  has a steeper minimum
    residual = abs(s_hermitian - s_hermitian['t2',::-1].runcopy(np.conj)).mean_all_but(['shift','t2'])
    residual.mean('t2')
    best_shift = residual.C.argmin('shift').item()

    s_FID = s['t2':(0,None)].C
    s_FID['t2',0] *= 0.5
    ph0 = s_FID['t2':0.0]
    ph0 /= abs(ph0)
    s_FID /= ph0
    if fl:
        fl.next('FID data')
        fl.image(s_FID.C.mean_all_but(['shift','t2']))
    s_FID.ft('t2')
    if fl:
        fl.next('FT of FID data')
        fl.image(s_FID.C.mean_all_but(['shift','t2']))
    logging.info(strm(ndshape(s_FID)))
    sum_abs_real = abs(s_FID.real).sum('t2').mean_all_but(['shift'])
    sum_abs_imag = abs(s_FID.imag).sum('t2').mean_all_but(['shift'])
    best_abs_real = (sum_abs_real/sum_abs_imag).argmin('shift').item()
    if fl:
        fig, (ax1,ax2,ax3) = plt.subplots(3,1)
        fl.next('mirror tests', fig=fig)
        def real_imag_mirror(forplot, ax):
            l = fl.plot(forplot.real, alpha=0.5, ax=ax, label='real',
                    human_units=False)
            fl.plot(forplot.C.setaxis('t2', lambda x: -x), color=l[-1].get_color(),alpha=0.2,ax=ax,
                    human_units=False)
            l=fl.plot(forplot.imag,alpha=0.5,ax=ax)
            fl.plot(-forplot.imag.C.setaxis('t2',lambda x: -x),color=l[-1].get_color(),
                    alpha=0.2,ax=ax,
                    human_units=False)
            residual = forplot['t2':(-max_shift,max_shift)]
            if ndshape(residual)['t2'] % 2 == 0:
                residual = residual['t2',:-1]
            fl.plot(abs(residual-residual['t2',::-1].runcopy(np.conj))*20, ax=ax, label='residual x 20',
                    human_units=False)
            ax.legend()
        #{{{show test for best_shift
        forplot = s['shift':best_shift].C.mean_all_but(['shift','t2'])
        forplot /= (forplot['t2':0]/abs(forplot['t2':0]))
        real_imag_mirror(forplot, ax1)
        ax1.set_title('best_shift')
        #}}}
        #{{{show test for shift=0
        forplot = s['shift':best_abs_real].C.mean_all_but(['shift','t2'])
        forplot /= (forplot['t2':0]/abs(forplot['t2':0]))
        real_imag_mirror(forplot,ax3)
        ax3.set_title('best abs real')
        #}}}
        #{{{ show test for shift = 0
        forplot = s['shift':0].C.mean_all_but(['shift','t2'])
        forplot /= (forplot['t2':0]/abs(forplot['t2':0]))
        real_imag_mirror(forplot, ax2)
        ax2.set_title('shift=0')
        #}}}
        fl.next('cost functions')
        def rescale(forplot):
            retval = forplot.C
            retval -= retval.data.min()
            retval /= retval.data.max()
            return retval
        fl.twinx(orig=True)
        residual.name('hermitian test')
        fl.plot(residual,c='k',
                human_units=False)
        fl.twinx(orig=False,color='red')
        fl.plot([],c='k',
                human_units=False)
        data_for_peak *= max_val
        data_for_peak.name('absolute value')
        fl.plot(data_for_peak.C.rename('t2','shift').set_units('shift','s'),c='red',
                human_units=False)
        xlim(np.array(rel_shift)*max_shift/1e-3)
    return best_shift+peak_center,max_shift    
