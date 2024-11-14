"""Rescaling QESR spectra
=========================
Demonstration of the rescaling of QESR spectra based on the
parameters of the dataset as well as the diameter and 
calibration info in the users pyspecdata config file
(``~/.pyspecdata`` on Linux/Mac or ``~/_pyspecdata`` on
Windows).
To run this example you should have the following parameters
set in your config file:
    :220720 propfactor: 8.637e-5
    :220720 q: 4600
    :220720 diameter: 0.704
"""
import pyspecdata as psd
from pyspecProcScripts import QESR_scalefactor
import matplotlib.pyplot as plt

file_loc = "francklab_esr/alex"
fig, ax_list = plt.subplots(1, 2)
# QESR_scalefactor takes the propfactor and diameter from my config file and I
# want to use the values that I've designated with the label 220720
calib_info = "220720"
with psd.figlist_var() as fl:
    d = psd.find_file("220922_70mM_270deg_teflon.DSC", exp_type=file_loc)
    if "harmonic" in d.dimlabels:
        d = d["harmonic", 0]
    fl.next("QESR rescaling", fig=fig)
    fl.plot(d, ax=ax_list[0])
    ax_list[0].set_title("Raw ESR spectra")
    d /= QESR_scalefactor(
        d, calibration_name=calib_info, diameter_name=calib_info
    )
    fl.plot(d, ax=ax_list[1])
    ax_list[1].set_title("Rescaled ESR spectra")
