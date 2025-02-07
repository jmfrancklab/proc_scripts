from pyspecdata.datadir import pyspec_config  # piggyback on _pyspecdata
from pyspecdata.general_functions import strm
import logging
from pint import UnitRegistry
from numpy import sqrt

default_Q = float(pyspec_config.get_setting("default Q", default=4700))
# note that this will fail if you don't have "QESR Conversion" set in your pyspecdata file

dint_conversion = pyspec_config.get_setting("QESR Conversion")
if dint_conversion is not None:
    dint_conversion = float(dint_conversion)

ureg = UnitRegistry(
    system="mks",
    autoconvert_offset_to_baseunit=True,
    auto_reduce_dimensions=True,
)
Q_ = ureg.Quantity


def QESR_scalefactor(d):
    "determine the scaling factor for QESR"
    if dint_conversion is None:
        raise ValueError("you need to set the parameter 'QESR Conversion' in you pyspecdata config file")
    # {{{ determine the signal denominator from the parameters of interest
    G_R = Q_(*d.get_prop("Gain"))
    C_t = Q_(*d.get_prop("ConvTime"))
    power = Q_(*d.get_prop("Power"))
    B_m = Q_(*d.get_prop("ModAmp"))
    Q = Q_(default_Q, "dimensionless")  # hard set Q value
    n_B = Q_(1, "dimensionless")  # calculate this
    S = Q_(0.5, "dimensionless")
    c = Q_(
        1, "dimensionless"
    )  # the first fraction on pg 2-17 -- essentially the conversion factor
    signal_denom = G_R * C_t * sqrt(power) * B_m * n_B * S * (S + 1) * Q
    signal_denom = signal_denom.to(Q_("G") * sqrt(Q_("W")) * Q_("s"))
    # }}}
    logging.debug(
        strm(
            f"$G_R={G_R:~P}$\n",
            f"$C_t={C_t:~P}$\n",
            f"$power={power:~P}$\n",
            f"$B_m={B_m:~P}$\n",
            f"$Q={Q:~P} $\n",
            f"$n_B={n_B:~P} $\n",
            f"$S={S:~P} $\n",
            f"$c={c:~P} $\n",
            f"signal denom$={signal_denom:~P}$",
            f"doubleint conversion$={dint_conversion}$",
        )
    )
    # normally, we divide by signal_denom.magnitude and multiply by dint_conversion and divide by 1e-6
    return signal_denom.magnitude / dint_conversion * 1e-6
