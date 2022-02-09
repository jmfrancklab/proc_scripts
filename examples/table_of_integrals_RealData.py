"""Convert 2D  Real Data into Integrals with Errors
===================================================

Take a 2D dataset and convert it to a table of integrals with errors, utilizing
all the bells and whistles (frequency and time selection, alignment, etc.)

Demonstrate on a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise, and with fluctuating field
(normally distributed field variation).
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from numpy.random import normal, seed
import sympy as s
from collections import OrderedDict

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 6

init_logging(level="debug")
fl = fl_mod()
signal_pathway = {"ph1": 1}
t_max = 0.1

with figlist_var() as fl:
    for (
        filename,
        exp_type,
        nodename,
        postproc,
        indirect,
        clock_correction,
        label,
        f_range,
    ) in [
        (
            "211223_Ras_M67R1a_capProbe",
            "ODNP_NMR_comp/ODNP",
            "enhancement",
            "spincore_ODNP_v1",
            "power",
            False,
            "M67 Ras mutation",
            (-0.5e3, 0.5e3),
        )
    ]:
        fl.basename = "(%s)" % label
        s = find_file(
            filename,
            exp_type=exp_type,
            expno=nodename,
            postproc=postproc,
            lookup=lookup_table,
        )
        s_int, s = peak_intensities(
            s,
            signal_pathway=signal_pathway,
            searchstr=label,
            f_range=f_range,
            t_max=t_max,
            indirect=indirect,
            Ep_real=True,
            alias_slop=2,
            clock_correction=clock_correction,
            fl=fl,
        )
