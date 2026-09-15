"""
FID from Echo after Phasing and Timing Correction
=================================================

Demonstrate how we can automatically find the zeroth order phase and the
center of the echo and then slice, in order to get a properly phased FID.
Here we see this

This example provides a relatively routine example.
"""

import pyspecdata as psd
import pyspecProcScripts as pypcs
import matplotlib.pyplot as plt
from pyspecProcScripts.load_data import lookup_table

psd.init_logging(level="info")
plt.rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 1

signal_pathway = {"ph1": 0, "ph2": 1}
with psd.figlist_var(file_name="tempdata220922final.pdf") as fl:
    for nodename in [
        "FIR_noPower",
        "FIR_noPower_1",
        "FIR_noPower_2",
        "FIR_noPower_3",
        "FIR_noPower_4",
    ]:
        d = psd.find_file(
            "260820_TTPDI_ODNP_3.h5",
            exp_type="B27/ODNP",
            postproc="spincore_IR_v4",
            lookup=lookup_table,
            expno=nodename,
        )  # returns signal with t=0 set approximately correctly
        fl.basename = nodename
        d = pypcs.fid_from_echo(d, signal_pathway, fl=fl)
        fl.next("final phased spectrum")
        fl.image(d)
        # in the following, I assume the units are auto-scaled to kHz
        print("inh_bounds", d.get_prop("inh_bounds"))
        plt.axvline(x=d.get_prop("inh_bounds")[0] / 1e3, color="w", ls=":")
        plt.axvline(x=d.get_prop("inh_bounds")[1] / 1e3, color="w", ls=":")
