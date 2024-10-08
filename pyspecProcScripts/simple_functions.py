"First order functions for very simple (a few lines) data manipulation"
import numpy as np
import logging

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
            return np.concatenate(
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


def find_apparent_anal_freq(s):
    """A function to identify the position of analytic signal as acquired on the oscilloscope.
    Importantly this function takes into account the effects of aliasing in identifying the
    frequency of the resulting signal.

    Parameters
    ==========
    s: nddata
        data with a single (dominant) peak, where you want to identify the frequency.

    Returns
    =======
    s: nddata
        The original data that was initially fed to the function
    nu_a: float
        The apparent frequency of the signal
    isflipped: boolean
        If aliased from a negative frequency, this notes whether the phase
        of the final time domain signal will be flipped
    """
    carrier = s.get_prop("acq_params")["carrierFreq_MHz"] * 1e6
    dt = s["t"][1] - s["t"][0]
    # In analytic signal the negative frequencies were tossed.
    # So we first need to check if the carrier is in the SW of
    # the analytic signal or if the signal was aliased
    if carrier < 1 / dt:
        logging.debug(
            "You are in the clear and no aliasing took place!"
        )
        nu_a = carrier
        isflipped = False
    else:
        logging.info(
            "Aliasing occurred, but we can still find that frequency!"
        )
        SW = 2 / dt  # SW of data before made analytic
        #              - what scope sees
        n = np.floor(  # nearest integer multiple of
            #            sampling frequency.
            #            Measured from the left side of
            #            the shifted spectrum
            (carrier + SW / 2)
            / SW
        )
        nu_a = carrier - n * SW
        isflipped = False
        if nu_a < 0:
            # when stored, we threw out the negative frequencies
            # so that means that we captured an aliased copy of the
            # negative frequency
            nu_a = -carrier + n * SW
            # we need to make note of the fact that the phase of our final time domain signal
            # (after filtering and up-conversion) will be flipped
            isflipped = True
    return s, nu_a, isflipped
