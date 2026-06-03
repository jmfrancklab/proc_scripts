"""Align a 2D ESR data set
=======================

This example shows how to align ESR spectra that share one direct field axis
and vary along an indirect dimension such as temperature or time.  The example
uses ``example.DSC`` as a placeholder; replace it with a real ESR file that
contains a ``Temperature`` dimension.
"""

from matplotlib.pyplot import rcParams
from pyspecdata import find_file, figlist_var
from pyspecProcScripts import align_esr_2d, cumulant_rms

rcParams["image.aspect"] = "auto"

# sphinx_gallery_thumbnail_number = 3

# Load a 2D ESR file whose direct field dimension is $B_0$ and whose indirect
# dimension is Temperature.  This placeholder filename is intended to be
# swapped for a real file on the user's computer.
d = find_file("Temperature_vs_Time_250428.h5", exp_type="rs_proc_data", expno="all_spectra")

with figlist_var() as fl:
    fl.next("raw 2D ESR")
    fl.image(d.real.C.reorder(["Temperature", "$B_0$"]))

    fl.next("cumulant RMS before alignment")
    fl.plot(cumulant_rms(d, "Time"), label="before")

    align_esr_2d(d, "Time", fl=fl)

    fl.next("aligned 2D ESR")
    fl.image(d.C.reorder(["Temperature", "$B_0$"]))

    fl.next("cumulant RMS after alignment")
    fl.plot(cumulant_rms(d, "Time"), label="after")
