from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import QESR_scalefactor
file_loc = 'francklab_esr/alex'
fieldaxis = '$B_0$'
background = find_file('220922_water.DSC', exp_type = file_loc)['harmonic',0]
with figlist_var() as fl:
    d = find_file('220922_70mM_270deg_teflon.DSC', exp_type = file_loc)
    if 'harmonic' in d.dimlabels:
        d = d['harmonic',0]
    color = d.get_plot_color()
    fl.next("Raw QESR")
    fl.plot(d,'o',color = color, alpha = 0.5)
    fl.next("Rescaled")
    d /= QESR_scalefactor(d, calibration_name = "220720")
    fl.plot(d)
    fl.show()
