"""Quantify the Double Integral of an ESR spectra
=================================================
Calculate the concentration of spin label present in
a sample based on the double integral of an ESR spectra.

An ESR spectra should be collected for both the 
sample of interest as well as a "background" 
(typically water for aqueous samples). 
The spectra is rescaled according to it's acquisition parameters 
(e.g., Q value, diameter of capillaries, modulation amplitude,
gain, etc) and multiplied by a proportionality constant that 
is defined in your pyspecdata config file. 

In order to properly run this example, your config file
should have the following values under General:

::

    220720 propfactor = 8.637e-5
    220720 q = 4600
    qesr caps diameter = 0.704
"""

from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import QESR_scalefactor

# sphinx_gallery_thumbnail_number = 2

background_file = find_file("220804_water.DSC", exp_type="francklab_esr/Farhana")["harmonic", 0]

with figlist_var() as fl:
    c = QESR(
        "220804_rasI36_MTSL.DSC", #filename
        label="kRas I36", #label for legends
        exp_type="francklab_esr/Farhana", #location of file
        calibration_name="220720", 
        diameter_name="QESR caps", 
        background=background_file, #background used for background subtraction - loaded above
        pickle_file="TEMPOL_rerun_conc.pickle",
        fl=fl
    )
