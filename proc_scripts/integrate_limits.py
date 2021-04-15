from pyspecdata import *
import numpy as np
import logging


def integrate_limits(s, axis="t2", convwidth=50):
    signal_sign = s.C.sum(axis).run(np.real).run(np.sign)
    logging.debug(strm("the signal sign", signal_sign))
    frq_slice = (
        (s.real * signal_sign)
        .mean_all_but(axis)
        .convolve(axis, convwidth)
        .contiguous(lambda x: abs(x) > 0.5 * abs(x).data.max())[0]
    )
    return frq_slice
