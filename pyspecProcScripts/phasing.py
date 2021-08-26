"""This module includes routines for phasing NMR spectra."""
from pyspecdata import *
from matplotlib.patches import Ellipse
from scipy.optimize import minimize
from pylab import *
import numpy as np
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
    eigenValues, eigenVectors = linalg.eigh(cov_mat)
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
def hermitian_function_test(s, selection_range=(None,0.04), ini_delay = 0e6, 
        band_mask_no = 0, band_mask = False, fl=None):
    r"""determine the center of the echo

    Parameters
    ==========
    selection_range:  tuple of floats
        Slice that extends to about 2 times the time domain signal.
    ini_delay:       float
        initial delay
    band_mask_no:    int
        Helps determing number of points for rectangular mask
    band_mask:       boolean
        Masks out aliasing components around the middle of the 
        selection range.
        When true, uses rectangular shape for applied mask.
        When false, uses triangular shape for applied mask.
    """
    orig_dt = s.getaxis('t2')[1] - s.getaxis('t2')[0]
    s.ft('t2')
    s = s['t2':(-5e3/2,5e3/2)] #assumes signal is somewhat on resonance
    s.ift('t2',pad=2048*8)
    if fl is not None:
        fl.next('data with padding')
        s = s['t2':(ini_delay,None)]
        fl.plot(abs(s),'.')
    selection = s['t2':selection_range]
    N = ndshape(selection)['t2']
    mid_idx = N//2+N%2-1
    selection = selection['t2',0:2*mid_idx+1]
    dt = selection.get_ft_prop('t2','dt')
    #the shifts themselves run from 0 to mid_idx -- the echo-centers
    #these correspond to are different. Also, we're not really
    #worried about the sub-integer shift here, because we can just
    #use sinc interpolation before we start -- see:
    #https://jmfrancklab.slack.com/archives/CLMMYDD98/p1623354066039100
    shifts =  nddata(dt*(r_[0:mid_idx]),'shift')
    shifts.set_units('shift','s')
    logger.info(strm("Length of shifts dimension:", ndshape(shifts)['shift']))
    residual = selection.C
    residual.ft('t2')
    residual *= np.exp(-1j*2*pi*shifts*residual.fromaxis('t2'))
    residual.ift('t2')
    logger.info(strm("Length of t2 dimension:", ndshape(selection)['t2']))
    assert ndshape(selection)['t2'] % 2 == 1, "t2 dimension *must* be odd, please check what went wrong."
    #{{{phase correct and weigh 'residual' by 1/(signal amplitude)
    center_point = residual['t2',mid_idx]
    residual /= center_point
    #}}}
    if fl is not None:
        fl.next('shifted and phased')
        if 'power' in s.dimlabels:
            fl.image(residual['power',-4])
        else:
            fl.image(residual)
    residual = abs(residual - residual['t2',::-1].runcopy(np.conj)).mean_all_but(['shift','t2'])
    if fl is not None:
        fl.next('Residual 2D')
        fl.image(residual)
    if band_mask:
        n_mask_pts = int(band_mask_no/dt)
        if n_mask_pts % 2:
            n_mask_pts -= 1
        residual['t2',mid_idx+n_mask_pts:] = 0
        residual['t2',:mid_idx-n_mask_pts] = 0
        title_str = 'rectangular mask'
    else:
        #Here we would do the calculation outlined for the triangle mask,
        #but with only two new variables -- A and B
        A = r_[-mid_idx:mid_idx+1]
        A = -abs(A)
        A += mid_idx
        B = residual.fromaxis('shift').data/dt
        A = A.reshape(-1,1)
        B = B.reshape(1,-1)
        mask = np.greater_equal(A,B)
        mask = np.multiply(mask,1)
        mask = mask.astype(float)
        norm = np.floor(mid_idx-B)
        norm = norm.astype(float)
        mask /= ((norm)*2)
        mask = nddata(mask,['t2','shift'])
        mask.setaxis('t2',residual.getaxis('t2'))
        mask.setaxis('shift',residual.getaxis('shift'))
        if fl is not None:
            fl.next('mask')
            fl.image(mask)
        residual *= mask
        title_str = 'triangular mask'
    if fl is not None:
        fl.next('masked residual')
        fl.image(residual)
    residual.rename('shift','center').setaxis('center',
            dt*(mid_idx-r_[0:mid_idx])+ini_delay)
    residual = residual['center',::-1]
    if fl is not None:
        fl.next('masked residual -- relabeled')
        fl.image(residual)
    residual.mean('t2')
    residual.set_units('center','s')
    best_shift = residual['center':(500e-5,None)].C.argmin('center').item()
    #slices first few points out as usually the artifacts give a minimum
    if fl is not None:
        axvline(x=best_shift*1e6, c='white', linestyle=":")
        fl.next('cost function %s - freq filter'%title_str)
        residual.name('cost function')
        fl.plot(residual, color='k', alpha = 0.5, human_units=False)
        fl.twinx(orig=False, color='red')
        selection.name('absolute value')
        fl.plot(abs(selection).mean_all_but('t2').rename('t2','center').set_units('center',
            's'), color='red',alpha=0.5, human_units=False)
        axvline(x=best_shift, c='k', linestyle=":")
        for dwell_int in r_[0:5]:
            axvline(x=best_shift-(orig_dt*dwell_int),alpha=0.4,c='k',linestyle=":")
        fl.twinx(orig=True)
        ylim(0,residual['center':(best_shift-4e-3,best_shift)].data.max())
    return best_shift    


