import pyspecdata as psd
from pyspecProcScripts import QESR_scalefactor
import matplotlib.pyplot as plt

file_loc = "francklab_esr/alex"
fig, ax_list = plt.subplots(1, 2)
with psd.figlist_var() as fl:
    d = psd.find_file("220922_70mM_270deg_teflon.DSC", exp_type=file_loc)
    if "harmonic" in d.dimlabels:
        d = d["harmonic", 0]
    fl.next("QESR rescaling", fig=fig)
    fl.plot(d, ax=ax_list[0])
    ax_list[0].set_title("Raw ESR spectra")
    d /= QESR_scalefactor(d, calibration_name="220720", diameter_name="220720")
    fl.plot(d, ax=ax_list[1])
    ax_list[1].set_title("Rescaled ESR spectra")
