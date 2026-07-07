"""
Phasing and Timing Correction Using a Varied Tau Experiment
===========================================================

Take real data with varying echo times, 
and demonstrate how we can automatically find the zeroth order phase and the
center of the echo in order to get data that's purely real in the frequency
domain.

Here, we specifically check to see whether or not the offset between the
programmed τ and the center of the echo
(as found by :func:`hermitian_function_test`)
is consistent.
"""
from pyspecdata import *
from pyspecProcScripts import *
from pylab import *
import sympy as s
from collections import OrderedDict

init_logging(level="debug")

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 1
t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
f_range = (-400, 400)
filename = "201113_TEMPOL_capillary_probe_var_tau_1"
signal_pathway = {"ph1": 1, "ph2": 0}

# TODO ☐: you are hacking the postprocessing!!! you can't just do that! You really should be able to use the default/existing postprocessing, but if that doesn't work, you need to carefully explain why this is needed
def process_var_tau_data(data):
    if "ph1" not in data.dimlabels:
        data.chunk("t", ["ph2", "ph1", "t2"], [2, 4, -1])
        data.setaxis("ph2", r_[0, 2] / 4)
        data.setaxis("ph1", r_[0:4] / 4)
    data.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    data.set_units("t2", "s")
    data *= data.get_prop("acq_params").get("nScans", 1)
    data.squeeze()
    data *= 2e-6 / 1.11e4
    data.set_units("V")
    data.ft("t2", shift=True).ft(["ph1", "ph2"])
    data.reorder(["ph1", "ph2", "tau"])
    return data


with figlist_var() as fl:
    for nodename, file_location, label in [
        (
            "var_tau",
            "ODNP_NMR_comp/test_equipment/var_tau",
            "tau is 1 ms",
        ),
    ]:
        data = find_file(
            filename,
            exp_type=file_location,
            expno=nodename,
        )
        # TODO ☐: the following is absolutely forbidden, and not within the scope of allowed changes!  You are adding a line to manually do what the postprocessing should do.
        data = process_var_tau_data(data)
        data = data["tau", :-7]
        tau_list = list(data.getaxis("tau"))
        data.reorder(["ph1", "ph2", "tau", "t2"])
        data = data["t2":f_range]
        mytable = []
        mytable.append(
            ["programmed tau / ms", "estimated tau / ms", "difference / ms"]
        )
        for j in range(len(tau_list)):
            tablerow = []
            alias_slop = 3
            programmed_tau = tau_list[j]
            tablerow.append(programmed_tau / 1e-3)
            this_data = data["tau", j]
            this_data.ift("t2")
            # TODO ☐: this is the same as what was here before except (1) it's MUCH more verbose (2) you unnecessarily create an intermediate variable and (3) you unset the units.  The first two are straight undesired.  It's not clear what's motivating the third
            for_shift = select_pathway(this_data, signal_pathway)
            for_shift.set_units(None)
            fl.basename = "%0.1f ms" % (programmed_tau / 1e-3)
            best_shift = hermitian_function_test(
                for_shift,
                aliasing_slop=alias_slop,
            )
            tablerow.append(best_shift / 1e-3)
            diff = abs(best_shift - programmed_tau)
            tablerow.append(diff / 1e-3)
            mytable.append(tablerow)

        def tabulate(mytable):
            print(" ".join(mytable[0]))
            strlens = [len(j) for j in mytable[0]]
            print(" ".join("-" * j for j in strlens))
            formatstr = " ".join(f"%{str(j)}.2f" for j in strlens)
            for j in mytable[1:]:
                print(formatstr % tuple(j))

        tabulate(mytable)
