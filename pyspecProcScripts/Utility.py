"These are simple utility functions"
from numpy import log10, ndarray, array


def dBm2power(dBm):
    "convert from dBm to W"
    if type(dBm) is not ndarray:
        dBm = array(dBm)
    return 10 ** (0.1 * dBm - 3)


def power2dBm(power):
    "convert from dBm to W"
    if type(power) is not ndarray:
        power = array(power)
    return log10(power) * 10 + 30


def Vpp2power(Vpp, Z=50.0):
    "convert peak to peak voltage to power"
    if type(Vpp) is not ndarray:
        Vpp = array(Vpp)
    Vrms = Vpp / 2 / sqrt(2)
    return Vrms**2 / Z
