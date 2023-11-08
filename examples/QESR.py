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
fieldaxis = '$B_0$'
with figlist_var() as fl:
    for filenum, (thisfile, that_exp_type, thislabel, pushout, calibration, diameter,
            background) in enumerate(
            [
              ('220804_rasI36_MTSL.DSC','francklab_esr/Farhana','kRas I36',0.4, '220720','QESR caps',
                  find_file('220804_water.DSC',exp_type='francklab_esr/Farhana')['harmonic',0]),
             ]):
        d = find_file(thisfile,exp_type = that_exp_type)
        d -= d[fieldaxis, -100:].data.mean()
        if "harmonic" in d.dimlabels:
            d = d['harmonic',0]
        color = d.get_plot_color()
        fl.next("Raw QESR")
        fl.plot(d,color=color, label=thislabel,alpha = 0.5)
        fl.next("Rescaled using the proportionality constant")
        d /= QESR_scalefactor(d, calibration_name = calibration,
                diameter_name = diameter)
        fl.plot(d,color=color,label=thislabel)
        c = protein_QESR(thisfile, label=thislabel,pushout=pushout,
                exp_type = that_exp_type,
                color=color,
                calibration_name = calibration,
                diameter_name = diameter,
                background =  background,
                which_plot = 'compare',
                pickle_file = "TEMPOL_rerun_conc.pickle",
                fl = fl)
