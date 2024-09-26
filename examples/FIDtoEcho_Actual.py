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

signal_pathway = {"ph1": 1}
with psd.figlist_var(file_name="tempdata220922final.pdf") as fl:
    for nodename in [
        "enhancement_10C",
        "enhancement_10C_repeat",
        "enhancement_15C",
        "enhancement_21C",
        "enhancement_5C",
    ]:
        d = psd.find_file(
            "211103_TEMPOL_269uM_HeatExch.h5",
            exp_type="ODNP_NMR_comp/ODNP",
            postproc="spincore_ODNP_v1",
            lookup=lookup_table,
            expno=nodename,
        )  # returns signal with t=0 set approximately correctly
        fl.basename = nodename
        d = pypcs.fid_from_echo(d, signal_pathway, fl=fl)
        fl.next("final phased spectrum")
        fl.image(d)
        # in the following, I assume the units are auto-scaled to kHz
        print("peakrange", d.get_prop("peakrange"))
        plt.axvline(x=d.get_prop("peakrange")[0] / 1e3, color="w", ls=":")
        plt.axvline(x=d.get_prop("peakrange")[1] / 1e3, color="w", ls=":")
