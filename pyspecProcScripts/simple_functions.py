"First order functions for very simple (a few lines) data manipulation"
import numpy as np
#from .phasing import zeroth_order_ph


class logobj(object):
    def __init__(self, array_len=1000):  # just the size of the buffer
        self.log_list = []
        # {{{ this is a structured array
        self.log_dtype = np.dtype(
            [("time", "f8"), ("Rx", "f8"), ("power", "f8"), ("cmd", "i8")]
        )
        self.log_array = np.empty(array_len, dtype=self.log_dtype)
        self.log_dict = {
            0: ""
        }  # use hash to convert commands to a number, and this to look up the meaning of the hashes
        # }}}
        self.currently_logging = False
        self.log_pos = 0
        self.array_len = array_len
        return

    @classmethod
    def from_group(cls, h5group):
        """initialize a new log object with data loaded from the h5py group
        h5group (factory method"""
        thislog = cls()
        thislog.__setstate__(h5group)
        return thislog

    @property
    def total_log(self):
        "the log is stored internally as a list of arrays -- here return a single array for the whole log"
        if hasattr(self, "_totallog"):
            return self._totallog
        else:
            return concatenate(
                self.log_list + [self.log_array[: self.log_pos]]
            )

    @total_log.setter
    def total_log(self, result):
        self._totallog = result

    def __setstate__(self, inputdict):
        in_hdf = False
        if "dictkeys" in inputdict.keys():
            self.log_dict = dict(
                zip(inputdict["dictkeys"], inputdict["dictvalues"])
            )
        elif "dictkeys" in inputdict.attrs.keys():
            # allows setstate from hdf5 node
            self.log_dict = dict(
                zip(inputdict.attrs["dictkeys"], inputdict.attrs["dictvalues"])
            )
            in_hdf = True
        else:
            raise IOError("I can't find dictkeys!")
        if in_hdf:
            self.total_log = inputdict["array"][
                :
            ]  # makes accessible after hdf is closed (forces into memory)
        else:
            self.total_log = inputdict["array"]


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


def determine_sign(s, signal_pathway, indirect, direct="t2", fl=None):
    """Determines the sign of the signal based on the difference between the 
    signal with the zeroth order phase correction applied to the entire data set vs 
    applying the zeroth order phase correction to each individual indirect indice.
    The sign can be used, so that when multiplied by the signal the data has a uniform 
    sign along the indirect axis.

    Parameters
    ==========
    s: nddata
        data with a single (dominant) peak, where you want to return the sign
        of the integral over all the data.
    signal_pathway: dict
        dictionary containing the coherence transfer pathway where the signal
        resides.
    indirect: str
        Name of the indirect axis along which the sign is being determined
    direct: str (default "t2")
        Name of the direct dimension

    Returns
    =======
    mysign: nddata
        A dataset with all +1 or -1 (giving the sign of the original signal).
        Does *not* include the `direct` dimension
    """
    assert s.get_ft_prop(
        direct
    ), "this only works on data that has been FT'd along the direct dimension"
    d = select_pathway(s,signal_pathway)
    d_integral = d.C.real.integrate(direct)
    ind_phase_ind = select_pathway(d.C,signal_pathway)
    for j in range(len(s.getaxis(indirect))):
        ind_phase_ind[indirect,j] /= zeroth_order_ph(ind_phase_ind[indirect,j])
    ind_phase_ind.run(np.real).integrate(direct)
    mysign = (ind_phase_int/d_integral).angle / np.pi
    mysign = np.exp(1j*np.pi*mysign.run(np.round))
    if fl is not None:
        fl.next("check sign")
        fl.image(d * mysign)
    return mysign
