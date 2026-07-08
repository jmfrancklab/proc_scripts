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
d = find_file("Temperature_vs_Time_250428.h5", exp_type="rs_proc_data", expno="all_spectra").real
d.reorder("$B_0$", first=False)

with figlist_var() as fl:
    fl.next("raw 2D ESR")
    fl.image(d)

    fl.next("difference from average, raw")
    fl.image(d['Time':(0,400)] - d['Time':(0,400)].C.mean('Time'))

    fl.next("difference from average, raw, plot a few")
    fl.plot((d['Time':(0,400)] - d['Time':(0,400)].C.mean('Time'))['Time',::100])

    fl.next("cumulant RMS", legend=True)
    fl.plot(cumulant_rms(d, "Time"), label="before alignment")

    align_esr_2d(d, "Time", fl=fl)

    fl.next("aligned 2D ESR")
    fl.image(d.reorder("$B_0$", first=False))

    fl.next("cumulant RMS")
    fl.plot(cumulant_rms(d, "Time"), label="after alignment")

    fl.next("difference from average, after align")
    fl.image(d['Time':(0,400)] - d['Time':(0,400)].C.mean('Time'))
