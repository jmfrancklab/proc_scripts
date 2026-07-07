"""
Echoes with varied tau lengths
==============================

Processes data which has varying lengths of tau. 
Demonstrates how to load a h5 file. 
"""
from pylab import *
from pyspecdata import *
import h5py as h5

rcParams["image.aspect"] = "auto"
# sphinx_gallery_thumbnail_number = 3

# TODO ☐: here you have the same issue where you're trying to hack the postproc.  This is not acceptable!!  If there is an issue with the existing postproc, discuss in the chat.
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
    for filename, expno, exp_type, frequency, f_range in [
        (
            "201209_Ni_sol_probe_var_tau_",
            "var_tau",
            "test_equipment/var_tau",
            14.89e6,
            (-13.5e3, 0),
        )
    ]:
        fl.basename = filename
        logger.info(strm("analyzing", filename))
        # {{{ JF wanted to see what dataset is called
        fullname = search_filename(filename, exp_type=exp_type, unique=True)
        with h5.File(fullname, "r") as fp:
            logger.info(strm(fp.keys()))
        # }}}
        # TODO ☐: again, all the following changes are undesired -- discuss in chat to thread the needle
        d = find_file(
            filename,
            exp_type=exp_type,
            expno=expno,
        )
        d = process_var_tau_data(d)
        d = d["t2":f_range]
        d = d["ph1", +1]["ph2", -2]
        d.ift("t2")
        fl.next("echoes")
        fl.plot(d.real, alpha=0.2, linewidth=0.5)
        fl.plot(abs(d), alpha=0.5, linewidth=1)
        NV = (
            250e-6 * 55.4 * 2 * N_A
        )  # 400 μL, 55.4 M water molecs, 2 spins/molec
        nu0 = frequency
        LambdaNMR = 1.55e-4  # 1 G/√W
        I = 0.5
        Vsignal = (
            LambdaNMR
            * NV
            * (gammabar_H * 2 * pi)
            * I
            * (I + 1)
            * (hbar * 2 * pi * nu0) ** 2
            * sqrt(50)
        )
        Vsignal /= 3 * k_B * (273 + 20)
        axhline(y=Vsignal, alpha=0.2)
        logger.info(strm("Vsignal expected", Vsignal))
        fl.plot(Vsignal, "x", label="theoretical signal at 0")
