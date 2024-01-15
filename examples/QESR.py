"""Quantify the Double Integral of an ESR spectra
=================================================
A ESR is collected for both the sample of interest as well
as a "background" (typically water for aqueous samples). 
The spectra is rescaled by it's acquisition parameters 
(e.g., Q value, diameter of capillaries, modulation amplitude,
gain, etc) and multiplied by a proportionality constant that 
is defined in your pyspecdata config file. In order to 
properly run this example, your config file should have 
the following values under General:
    220720 propfactor = 8.637e-5
    220720 q = 4600
    qesr caps diameter = 0.704
"""

from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import QESR_scalefactor

# sphinx_gallery_thumbnail_number = 2

fieldaxis = "$B_0$"
background_file = find_file("220804_water.DSC", exp_type="francklab_esr/Farhana")["harmonic", 0]

with figlist_var() as fl:
    c = QESR(
        "220804_rasI36_MTSL.DSC",
        label="kRas I36",
        pushout=0.4,
        exp_type="francklab_esr/Farhana",
        calibration_name="220720",
        diameter_name="QESR caps",
        background=background_file,
        which_plot="compare",
        pickle_file="TEMPOL_rerun_conc.pickle",
        fl=fl,
    )
