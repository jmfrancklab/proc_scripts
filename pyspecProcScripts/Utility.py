"These are simple utility functions"

import numpy as np


def dBm2power(dBm):
    "convert from dBm to W"
    if type(dBm) is not np.ndarray:
        dBm = np.array(dBm)
    return 10 ** (0.1 * dBm - 3)


def power2dBm(power):
    "convert from dBm to W"
    if type(power) is not np.ndarray:
        power = np.array(power)
    return np.log10(power) * 10 + 30


def Vpp2power(Vpp, Z=50.0):
    "convert peak to peak voltage to power"
    if type(Vpp) is not np.ndarray:
        Vpp = np.array(Vpp)
    Vrms = Vpp / 2 / np.sqrt(2)
    return Vrms**2 / Z
