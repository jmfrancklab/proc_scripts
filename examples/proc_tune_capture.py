"""
Processing the Captured Tuning Curve
====================================

Takes the npz file of the captured tuning curve at different
zoom levels and plots them on the same plot, allowing us to 
look at any drift or discrepancies of the tuning curve.
"""
from pyspecdata import *
import matplotlib.pyplot as plt

fl = figlist_var()
for filename, thislabel in [
    ("220808_150uM_TEMPOL.npz", "150 uM"),
    ("220114_3mM_TEMPOL_3b.npz", "std 3 mM") #this is our control so don't change! You want your tuning curve to match this one
]:
    # {{{Load in npz file
    thisfile = search_filename(
        filename,
        exp_type="francklab_esr/alex",  # adjust the exp_type according to your personal folder
        unique=True,
    )
    data = np.load(thisfile)
    # }}}
    # {{{Go over range of zooms and plot against control
    for j in range(3):
        if j == 0:
            pass
        else:
            nd_data = {}
            zoom_data = data["zoom%d" % j].squeeze()
            zoom_data_nd = nddata(zoom_data[0], "frequency")
            zoom_data_nd.setaxis("frequency", zoom_data[1])
            nd_data["zoom%d" % j] = zoom_data_nd
            shift_val = zoom_data_nd.argmin()
            zoom_data_nd.setaxis("frequency", lambda x: x - shift_val["frequency"])
            fl.next("Tuning curve comparison-zoom%d" % j)
            fl.plot(zoom_data_nd, label=thislabel, alpha=0.5)
    # }}}
fl.show()
