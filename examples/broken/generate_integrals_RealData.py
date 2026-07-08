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
        # {{{ local replacement for the removed generate_integrals helper
        #     Git history: d9a49fe1 deleted
        #     pyspecProcScripts.generate_integrals because
        #     rough_table_of_integrals had taken over the full
        #     correction workflow.
        #     This old real-data example still demonstrates the final
        #     integration-with-error step, so the setup that the old
        #     helper bundled is kept locally before calling
        #     frequency_domain_integral below.
        if clock_correction:
            s = clock_correct(s, indirect=indirect, fl=fl)
        phase_dims = [j for j in s.dimlabels if j.startswith("ph")]
        if any(
            s.get_ft_prop(j) and not s.get_ft_prop(j, "unitary")
            for j in phase_dims
        ):
            # {{{ convert legacy non-unitary phase transforms
            #     Legacy spincore_ODNP_v1 used a non-unitary phase-cycle
            #     FT.  calc_masked_variance averages noise across
            #     phase-cycle axes and requires unitary phase transforms
            #     for that propagation.
            s.ift(phase_dims)
            for j in phase_dims:
                s.set_ft_prop(j, "unitary", True)
            s.ft(phase_dims)
            # }}}
        zero_pathway = {j: 0 for j in signal_pathway}
        excluded_pathways = [signal_pathway]
        if zero_pathway != signal_pathway:
            excluded_pathways.append(zero_pathway)
        s = s["t2":f_range]
        # }}}
        # See the d9a49fe1 note above: the old wrapper is gone, and this
        # is the lower-level helper that performs the integration/error
        # step.
        s_int = frequency_domain_integral(
            s,
            signal_pathway=signal_pathway,
            excluded_pathways=excluded_pathways,
            indirect=indirect,
            fl=fl,
        )
        fl.next(f"{label} integrals")
        fl.plot(s_int, ".", capsize=6)
