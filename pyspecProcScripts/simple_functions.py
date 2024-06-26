"First order functions for very simple (a few lines) data manipulation"
import numpy as np


def select_pathway(*args, **kwargs):
    r"""select a particular CT pathway from the signal `s`

    Arguments are *either* ``pathway`` -- a dict of key/value pairs indicating
    the pathway **or** the same set of key/value pairs, just passed as a dict.

    Parameters
    ==========
    s: nddata
        the data whose coherence pathway you would like to select

    pathway: dict
        keys are the names of the coherence transfer dimensions (conj. of phase
        cycling dimensions) and values are the pathway you want to select
    """
    if len(args) == 2 and len(kwargs) == 0:
        s, pathway = args
    elif len(args) == 1 and len(kwargs) > 0 and len(kwargs) % 2 == 0:
        s = args[0]
        pathway = kwargs
    else:
        raise ValueError("your arguments don't make any sense!!")
    retval = s
    for k, v in pathway.items():
        retval = retval[k, v]
    return retval


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
    assert s.get_ft_prop(
        direct
    ), "this only works on data that has been FT'd along the direct dimension"
    if fl is not None:
        fl.push_marker()
        fl.next("selected pathway")
        if "vd" in s.dimlabels:
            fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"))
        else:
            fl.image(s)
    data_sgn = s.C.sum(direct)
    data_sgn /= data_sgn.max().item()
    data_sgn.run(np.real).run(lambda x: np.sign(x))
    if fl is not None:
        fl.next("check sign")
        if "vd" in s.dimlabels:
            fl.image(
                s.C.setaxis("vd", "#").set_units("vd", "scan #") * data_sgn
            )
        else:
            fl.image(s * data_sgn)
        fl.pop_marker()
    return data_sgn
