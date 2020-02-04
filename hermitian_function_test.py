from pyspecdata import *
def zeroth_order_ph(d, plot_name=None):
    r'''determine the covariance of the datapoints
    in complex plane, and use to phase the
    zeroth-order even if the data is both negative
    and positive'''
    eigenValues, eigenVectors = eig(cov(c_[
        d.data.real,
        d.data.imag].T
        ))
    # next 3 lines from stackexchange -- sort by
    # eigenvalue
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    # determine the phase angle from direction of the
    # largest principle axis
    ph0 = arctan2(eigenVectors[1,0],eigenVectors[0,0])
    if plot_name:
        eigenVectors *= (eigenValues.reshape(-1,2)*ones((2,1)))/eigenValues.max()*abs(d.data).max()
        d_forplot = d.C
        fl.next(plot_name)
        fl.plot(
                d_forplot.data.real,
                d_forplot.data.imag,
                '.',
                alpha=0.25,
                label='before'
                )
        d_forplot /= exp(1j*ph0)
        fl.plot(
                d_forplot.data.real,
                d_forplot.data.imag,
                '.',
                alpha=0.25,
                label='after'
                )
        fl.plot(0,0,'ko')
        fl.plot(eigenVectors[0,0],eigenVectors[1,0],'o',
                label='first evec')
        fl.plot(eigenVectors[0,1],eigenVectors[1,1],'o')
    return exp(1j*ph0)
def hermitian_function_test(s, down_from_max=0.5):
    r"""determine the center of the echo"""
    # {{{ determine where the "peak" of the echo is,
    # and use it to determine the max
    # shift
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
    s_foropt = s.C
    # {{{ make sure there's an odd number of points
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
    n_points = ndshape(s_foropt)['t2']
    if n_points % 2 == 0:
        s_foropt = s_foropt['t2',:-1]
        n_points -= 1
    # }}}
    # {{{ construct test arrays for T2 decay and shift
    shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
    t2_decay = exp(-s_foropt.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
    # }}}
    # {{{ time shift and correct for T2 decay
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*
            s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    s_foropt /= t2_decay
    # }}}
    # the problem here is that we need to do this separately for each choice
    # of decay and shift, but allow zeroth_order_ph to operate over any
    # additional dimensions, so the following doesn't work correctly right now
    for j in range(ndshape(s_foropt)['R2']):
        for k in range(ndshape(s_foropt)['shift']):
            # {{{ phase the center point to predominantly 0 or pi
            center_point = s_foropt['R2',j]['shift',k]['t2',n_points//2+1]
            if len(center_point.dimlabels) > 0:
                assert len(center_point.dimlabels) == 1, repr(ndshape(center_point))+" has too many dimensions"
                ph0 = zeroth_order_ph(center_point)
                s_foropt['R2',j]['shift',k] /= ph0
            else:
                ph0 = center_point/abs(center_point)
                s_foropt['R2',j]['shift',k] /= ph0
            # }}}
    residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
    residual.reorder('shift')
    minpoint = residual.argmin()
    best_shift = minpoint['shift']
    best_R2 = minpoint['R2']
    return residual,best_shift+peak_center,best_R2
