from pyspecdata import *
def hermitian_function_test(s, down_from_max=0.5):
    # {{{ determine where the "peak" of the echo is,
    # center on it, and use it to determine the max
    # shift
    data_for_peak = abs(s).mean_all_but(['t2'])
    max_val = data_for_peak.data.max()
    pairs = data_for_peak.contiguous(lambda x: abs(x) >
            max_val*down_from_max)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    s.setaxis('t2',lambda x: x-peak_location.mean())
    s.register_axis({'t2':0})
    max_shift = diff(peak_location).item()/2
    # }}}
    # {{{ construct test arrays for T2 decay and shift
    shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
    t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
    # }}}
    s_foropt = s.C
    # {{{ time shift and correct for T2 decay
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*
            s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    s_foropt /= t2_decay
    # }}}
    # {{{ make sure there's and odd number of poings
    # and set phase of center point to 0
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
    n_points = ndshape(s_foropt)['t2']
    if n_points % 2 == 0:
        s_foropt = s_foropt['t2',:-1]
        n_points -= 1
    center_point = s_foropt['t2',n_points//2+1]
    s_foropt /= center_point/abs(center_point)
    # }}}
    residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
    residual.reorder('shift')
    minpoint = residual.argmin()
    best_shift = minpoint['shift']
    best_R2 = minpoint['R2']
    return residual,best_shift, best_R2
