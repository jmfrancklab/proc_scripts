from pyspecdata.datadir import pyspec_config # piggyback on _pyspecdata
from pint import UnitRegistry
from numpy import sqrt
class calib_info (object):
    "class for caching calibration info"
    def __init__(self):
        self.current_calib_name = None
        self.current_diam_name = None
    def use_calibration(self, calibration_name):
        if calibration_name != self.current_calib_name:
            assert type(calibration_name) is str, "the calibration name must be specified as a string!"
            assert len(calibration_name)>0, "you MUST use a calibration name and set the values `[calibration_name] q` and `[calibration_name] conversion` in your _pyspecdata file!"
            self.current_calib_name = calibration_name
            # note that this will fail if you don't have "[calibration_name] conversion" and "[calibration_name] q" set in your pyspecdata file
            self.default_Q = float(pyspec_config.get_setting(f"{calibration_name} Q"))
            self.dint_conversion = float(pyspec_config.get_setting(f"{calibration_name} Conversion"))
    def use_diameter(self, diameter_name):
        if diameter_name != self.current_diam_name:
            assert type(diameter_name) is str, "the diameter name must be specified as a string!"
            assert len(diameter_name)>0, "you MUST use a diameter name and set the values `[diameter_name] q` and `[diameter_name] conversion` in your _pyspecdata file!"
            self.current_diam_name = diameter_name
            self.d = float(pyspec_config.get_setting(f"{diameter_name} diameter"))
calibcache = calib_info()
ureg = UnitRegistry(
    system="mks", autoconvert_offset_to_baseunit=True, auto_reduce_dimensions=True
)
Q_ = ureg.Quantity
def QESR_scalefactor(d, calibration_name=None, diameter_name=None):
    "determine the scaling factor for QESR"
    calibcache.use_calibration(calibration_name)
    calibcache.use_diameter(diameter_name)
    # {{{ determine the signal denominator from the parameters of interest
    G_R = Q_(*d.get_prop("Gain"))
    C_t = Q_(*d.get_prop("ConvTime"))
    power = Q_(*d.get_prop("Power"))
    B_m = Q_(*d.get_prop("ModAmp"))
    Q = Q_(calibcache.default_Q, "dimensionless")  # hard set Q value
    d = Q_(calibcache.d, "mm")  # hard set Q value
    n_B = Q_(1, "dimensionless")  # calculate this
    S = Q_(0.5, "dimensionless")
    c = Q_(
        1, "dimensionless"
    )  # the first fraction on pg 2-17 -- essentially the conversion factor
    signal_denom = G_R * C_t * sqrt(power) * B_m * n_B * S * (S + 1) * Q * d**2
    signal_denom = signal_denom.to(Q_("G") * sqrt(Q_("W")) * Q_("s") * Q_("m")**2)
    # }}}
    print(
            f"$G_R={G_R:~L}$\n",
            f"$C_t={C_t:~L}$\n",
            f"$power={power:~L}$\n",
            f"$B_m={B_m:~L}$\n",
            f"$Q={Q:~L} $\n",
            f"$n_B={n_B:~L} $\n",
            f"$S={S:~L} $\n",
            f"$c={c:~L} $\n",
            f"$d={d:~L} $\n",
            f"signal denom$={signal_denom:~L}$",
            f"doubleint conversion$={calibcache.dint_conversion}$",
            )
    # normally, we divide by signal_denom.magnitude and multiply by calibcache.dint_conversion and divide by 1e-6
    return signal_denom.magnitude/calibcache.dint_conversion*1e-6
