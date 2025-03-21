"""QESR of same sample
===================

Here, we sit and acquire the same sample with different acquisition
parameters, and make sure we get the same QESR.

We do this b/c the data is not stored by XEPR in a way that we would
expect.
"""
from matplotlib.pyplot import axvline, axhline, gca
from pint import UnitRegistry
from pyspecdata import *
from itertools import cycle
from pyspecProcScripts import *
from pyspecProcScripts import QESR_scalefactor, QESR

init_logging(level='debug')

colors = plt.rcParams[
    "axes.prop_cycle"
]()  # this is the default matplotlib cycler for line styles
fieldaxis = "$B_0$"
plot_rescaled = False
water_nddata = find_file("230511_water.DSC", exp_type="francklab_esr/romana")["harmonic", 0]
with figlist_var() as fl:
    for file_searchstring, thislabel, exp_type in [
        (
            "250314_20mM_TSO4_water_AG_Parameter.DSC",
            "AG's parameter",
            "francklab_esr/romana",
        ),
        (
            "250314_20mM_TSO4_water_convTime4ms.DSC",
            "convTime 4 ms",
            "francklab_esr/romana",
        ),
        (
            "250314_20mM_TSO4_water_convTime8ms.DSC",
            "convTime 8 ms",
            "francklab_esr/romana",
        ),
        (
            "250314_20mM_TSO4_water_convTime12ms.DSC",
            "convTime 12 ms",
            "francklab_esr/romana",
        ),
        (
            "250314_20mM_TSO4_water_convTime24ms.DSC",
            "convTime 24 ms",
            "francklab_esr/romana",
        ),
    ]:
        logging.debug(strm("now analyzing",file_searchstring))
        c = QESR(
            file_searchstring,
            label=thislabel,
            exp_type=exp_type,
            plot_derivative=True,
            fl=fl,
        )
        print(f"{c:~P}")
