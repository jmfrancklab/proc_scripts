"First order functions for very simple (a few lines) data manipulation"
import numpy as np
def select_pathway(s,pathway):
    r"""select a particular CT pathway from the signal `s`
    
    Parameters
    ==========
    pathway: dict
        keys are the names of the coherence transfer dimensions (conj. of phase
        cycling dimensions) and values are the pathway you want to select
    """
    retval = s
    for k, v in pathway.items():
        retval = retval[k,v]
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
    assert s.get_ft_prop(direct), "this only works on data that has been FT'd along the direct dimension"
    if fl is not None:
        fl.push_marker()
        fl.next('selected pathway')
        fl.image(s)
    data_sgn = s.C.sum(direct)
    data_sgn /= data_sgn.max().item()
    data_sgn.run(np.real).run(lambda x: np.sign(x))
    if fl is not None:
        fl.next('check sign')
        fl.image(s*data_sgn)
        fl.pop_marker()
    return data_sgn
