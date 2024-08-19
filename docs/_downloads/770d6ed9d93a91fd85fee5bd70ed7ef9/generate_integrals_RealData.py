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
import sympy as s

init_logging(level="debug")

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 3

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
        s_int, s = generate_integrals(
            s,
            signal_pathway=signal_pathway,
            searchstr=label,
            f_range=f_range,
            indirect=indirect,
            alias_slop=2,
            clock_correction=clock_correction,
            fl=fl,
        )
