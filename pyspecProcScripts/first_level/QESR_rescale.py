from pyspecdata.datadir import pyspec_config  # piggyback on _pyspecdata
from pint import UnitRegistry
from numpy import sqrt
from pyspecdata import strm
import logging

logger = logging.getLogger("pySpecProcScripts.first_level.QESR_rescale")


class calib_info(object):
    "class for caching calibration info"

    def __init__(self):
        self.current_calib_name = None
        self.current_diam_name = None

    def use_calibration(self, calibration_name):
        """if diameter_name is None, default to the calibration name given
        by default calibration"""
        if calibration_name == "default":
            raise ValueError(
                "you may not use a calibration named 'default' -- that is"
                " just too confusing!!!\nGive your calibration a meaningful"
                ' name, and then set "default calibration" which gives the'
                " name that is used for the propFactor and for the Q factor as"
                ' well as "default diameter" which gives the name that is used'
                " for the diameter"
            )
        calibration_name_setting = calibration_name
        if calibration_name is None:
            calibration_name = pyspec_config.get_setting("default calibration")
        if calibration_name != self.current_calib_name:
            assert (
                type(calibration_name) is str
            ), "the calibration name must be specified as a string!"
            assert len(calibration_name) > 0, (
                "you MUST use a calibration name and set the values"
                " `[calibration_name] q` and `[calibration_name] propFactor`"
                " in your _pyspecdata file!"
            )
            # note that this will fail if you don't have
            # "[calibration_name] propFactor" and "[calibration_name] q"
            # set in your pyspecdata file
            try:
                self.default_Q = float(
                    pyspec_config.get_setting(f"{calibration_name} Q")
                )
                self.dint_propFactor = float(
                    pyspec_config.get_setting(f"{calibration_name} propFactor")
                )
            except Exception:
                raise RuntimeError(
                    "You set calibration_name to"
                    f' "{calibration_name_setting}"\nSo, I expect a line in'
                    " the [General] block of your pyspecdata config file (in"
                    " your home directory) that sets the value of the"
                    " following"
                    f" variables:\n{calibration_name} q\n"
                    f"{calibration_name} propFactor\n\n(I"
                    " can't run without these)"
                )
            self.current_calib_name = calibration_name

    def use_diameter(self, diameter_name):
        """if diameter_name is None, default to the calibration name given
        by default calibration"""
        if diameter_name is None:
            diameter_name = pyspec_config.get_setting("default diameter")
        if diameter_name != self.current_diam_name:
            assert (
                type(diameter_name) is str
            ), "the diameter name must be specified as a string!"
            assert len(diameter_name) > 0, (
                "you MUST use a diameter name and set the values"
                " `[diameter_name] q` and `[diameter_name] propFactor` in your"
                " _pyspecdata file!"
            )
            try:
                self.d = float(
                    pyspec_config.get_setting(f"{diameter_name} diameter")
                )
            except Exception:
                raise RuntimeError(
                    "(note this the second in a similar pair of errors) I"
                    " expect a line in the [General] block of your pyspecdata"
                    " config file (in your home directory) that sets the"
                    " value of the following"
                    f" variables:\n{diameter_name} diameter"
                )
            self.current_diam_name = diameter_name


calibcache = calib_info()
ureg = UnitRegistry(
    system="mks",
    autoconvert_offset_to_baseunit=True,
    auto_reduce_dimensions=True,
)
Q_ = ureg.Quantity


def QESR_scalefactor(d, **kwargs):
    raise ValueError("""QESR_scalefactor no longer exists
    replace

    d /= QESR_scalefactor(d)

    with 

    d = QESR_apply_scalefactor(d)

    and set the keyword arguments (calibration_name, diameter_name) as
    properties of the data:

    d.set_prop("calibration_name","blabla")
    d = QESR_apply_scalefactor(d)

    We are doing this because QESR_apply_scalefactor makes sure the scalefactor
    isn't applied more than once, and the properties are also used by
    align_ESR.
                     """)


def QESR_apply_scalefactor(d):
    """Divide the ESR spectrum by this number so that the double integral
    should be equal to concentration in μM.

    We specifically calculate
    :math:`d^{2} G_{R}  C_{t}  \\sqrt{P}  B_{m}  Q  n_{B}  S  (S + 1) /
    c_{propfactor}`
    (the denominator of equation 2-17 in the E500
    manual with :math:`d^2` (diameter squared) added in,
    so that after dividing by this term,
    we get a concentration of spins
    rather than the total number of spins).

    The constant :math:`c_{propfactor}` is determined by the calibration name.
    Note that a larger :math:`c_{propfactor}` corresponds to a higher
    concentration
    for the same recorded spectrum.

    Parameters
    ==========
    d : nddata
        has these properties:

        calibration_name:   str
                            The key corresponding to the appropriate
                            proportionality constant
                            in your pyspecdata config file.
                            Typically this is one value that doesn't need
                            changing
        diameter_name:      str
                            The key corresponding to the diameter of the
                            capillary tube used in the ESR experiment. This
                            enables us to calculate a reliable concentration
                            regardless of which capillary was used.

    Returns
    =======
    Denominator of equation 2-17 divided by the proportionality constant.
    Also includes a factor of 1e-6 to convert Molar to micromolar.
    """
    if d.get_prop("has_been_calibrated"):
        raise ValueError("this spectrum has already been calibrated!!")
    else:
        d.set_prop("has_been_calibrated", True)
    calibcache.use_calibration(d.get_prop("calibration_name"))
    calibcache.use_diameter(d.get_prop("diameter_name"))
    # {{{ determine the signal denominator from the parameters of interest
    G_R = Q_(*d.get_prop("Gain"))
    C_t = Q_(*d.get_prop("ConvTime"))
    power = Q_(*d.get_prop("Power"))
    B_m = Q_(*d.get_prop("ModAmp"))
    Q = Q_(calibcache.default_Q, "dimensionless")  # hard set Q value
    diameter = Q_(calibcache.d, "mm")  # diameter
    n_B = Q_(1, "dimensionless")  # calculate this
    S = Q_(0.5, "dimensionless")
    c = Q_(
        1, "dimensionless"
    )  # the first fraction on eq 2-17 -- in bruker E500 manual
    c_propfactor = Q_(calibcache.dint_propFactor, "m**2")
    dint_conversion = (c_propfactor / diameter**2).to("").magnitude
    signal_denom = (
        G_R * Q_(4, "ms") * sqrt(power) * B_m * n_B * S * (S + 1) * Q
    )
    signal_denom = signal_denom.to(Q_("G") * sqrt(Q_("W")) * Q_("s"))
    # }}}
    logger.debug(
        strm(
            f"$G_R={G_R:~P}$\n"
            f"$C_t={C_t:~P}$\n"
            f"$power={power:~P}$\n"
            f"$B_m={B_m:~P}$\n"
            f"$Q={Q:~P} $\n"
            f"$n_B={n_B:~P} $\n"
            f"$S={S:~P} $\n"
            f"$c={c:~P} $\n"
            f"$d={diameter:~P} $\n"
            f"signal denom$={signal_denom:~P}$\n"
            f"doubleint propFactor$={calibcache.dint_propFactor}$\n"
            f"doubleint conversion$={dint_conversion}$\n"
        )
    )
    # normally, we divide by signal_denom.magnitude and multiply by
    # calibcache.dint_propFactor and divide by 1e-6
    d *= dint_conversion / signal_denom.magnitude / 1e-6
    return d
