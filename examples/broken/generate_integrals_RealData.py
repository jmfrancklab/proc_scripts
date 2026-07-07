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
        # TODO ☐: you have added an entire code block here without adding an explanation!  You need to explain!
        # TODO ☐: it seems quite likely that you need to research the code history in git, and explain what what changed and what has motivated your changes here (when I tell you to explain, I mean with inline comments in the relevant location!)
        if clock_correction:
            s = clock_correct(s, indirect=indirect, fl=fl)
        phase_dims = [j for j in s.dimlabels if j.startswith("ph")]
        if any(
            s.get_ft_prop(j) and not s.get_ft_prop(j, "unitary")
            for j in phase_dims
        ):
            s.ift(phase_dims)
            for j in phase_dims:
                s.set_ft_prop(j, "unitary", True)
            s.ft(phase_dims)
        zero_pathway = {j: 0 for j in signal_pathway}
        excluded_pathways = [signal_pathway]
        if zero_pathway != signal_pathway:
            excluded_pathways.append(zero_pathway)
        s = s["t2":f_range]
        # TODO ☐: you have changed the name of the function called here, without adding a comment to explain what's going on
        s_int = frequency_domain_integral(
            s,
            signal_pathway=signal_pathway,
            excluded_pathways=excluded_pathways,
            indirect=indirect,
            fl=fl,
        )
        fl.next(f"{label} integrals")
        fl.plot(s_int, ".", capsize=6)
