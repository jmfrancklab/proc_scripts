"""Used to preprocess data based on type of experiment performed. Returns all data FTed
into the frequency domain with the phase cycles also FTed (coherence domain). Data is
not sliced or altered in anyway.

Parameters
==========
self:   nddata or h5 file

Returns
=======
nddata that has been FTed and in coherence domain
"""
from pyspecdata import *
from .Utility import dBm2power
import os
from sympy import symbols
import logging
import numpy as np
import logging
from pylab import *
from .DCCT_func import DCCT


def proc_spincore_SE_v1(s, fl=None):
    s.ft("ph1")  # In order to match the amplitude of the
    #             coherence domain signal to its size in
    #             each transient we do not use unitary FT
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    return s


def proc_spincore_diffph_SE_v1(s, fl=None):
    s = proc_spincore_diffph_SE_v2(s, fl=fl)
    s *= s.get_prop("acq_params")["nScans"]
    return s


def proc_spincore_diffph_SE_v2(s, fl=None):
    r"""this one uses a phase cycle where the overall phase and 90-180
    phase difference are cycled in a nested way -- see the DCCT paper to
    understand this!"""
    s.ft(["ph2", "ph_diff"])  # if we have used cycles for the axis
    #                          coordinates, signal in the coherence
    #                          dimension will match the amplitude of signal
    #                          in a single transient if we do this
    # {{{ after the FT, these have a different meaning in terms of coherence
    #     pathways -- remember that when labeling, pySpecData will change
    #     the ph here to a δp
    s.rename("ph2", "ph_overall")  # overall change in coherence
    s.rename("ph_diff", "ph1")  # change during pulse 1
    # }}}
    s.set_prop("coherence_pathway", {"ph_overall": -1, "ph1": +1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s.squeeze()
    s.ft("t2", shift=True)
    s.reorder(["ph1", "ph_overall"])
    return s


def proc_Hahn_echoph(s, fl=None):
    logging.debug("loading pre-processing for Hahn_echoph")
    nPhaseSteps = 8
    SW_kHz = s.get_prop("acq_params")["SW_kHz"]
    nScans = s.get_prop("acq_params")["nScans"]
    s.reorder("t", first=True)
    s.chunk("t", ["ph2", "ph1", "t2"], [2, 4, -1])
    s.labels({"ph2": r_[0.0, 2.0] / 4, "ph1": r_[0.0, 1.0, 2.0, 3.0] / 4})
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.reorder(["ph2", "ph1"])
    s.setaxis("nScans", r_[0:nScans])
    s.reorder("t2", first=False)
    s.ft("t2", shift=True)
    if fl is not None:
        fl.next("raw data, chunked")
        fl.image(abs(s))
    s.ft(["ph1", "ph2"])
    if fl is not None:
        fl.next("coherence")
        fl.image(abs(s))
    return s


def proc_spincore_IR(s, fl=None):
    vd_axis = s.getaxis("vd")
    if "t" in s.dimlabels:
        s.chunk("t", ["ph2", "ph1", "t2"], [2, 2, -1])
    s.setaxis("ph1", r_[0, 2.0] / 4)
    s.setaxis("ph2", r_[0, 2.0] / 4)
    s.reorder(["ph1", "ph2"]).set_units("t2", "s")
    s.set_prop("coherence_pathway", {"ph1": 0, "ph2": +1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft(["ph1", "ph2"])
    if fl is not None:
        fl.next("raw data -- coherence channels")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"))
    s.ift("t2")
    if fl is not None:
        fl.next("time domain (all $\\Delta p$)")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"))
    s.ft("t2")
    if fl is not None:
        fl.next("frequency domain (all $\\Delta p$)")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"), black=False)
    return s


def proc_spincore_IR_v2(s, fl=None):
    vd_axis = s.getaxis("vd")
    if "t" in s.dimlabels:
        s.chunk("t", ["ph2", "ph1", "t2"], [4, 4, -1])
    s.setaxis("ph1", r_[0, 1, 2, 3.0] / 4)
    s.setaxis("ph2", r_[0, 1, 2, 3.0] / 4)
    s.set_prop("coherence_pathway", {"ph1": 0, "ph2": -1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.reorder(["ph1", "ph2"]).set_units("t2", "s")
    s.ft("t2", shift=True)
    s.ft(["ph1", "ph2"])
    if fl is not None:
        fl.next("raw data -- coherence channels")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"))
    s.ift("t2")
    if fl is not None:
        fl.next("time domain (all $\\Delta p$)")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"))
    s.ft("t2")
    if fl is not None:
        fl.next("frequency domain (all $\\Delta p$)")
        fl.image(s.C.setaxis("vd", "#").set_units("vd", "scan #"), black=False)
    return s


def proc_nutation(s, fl=None):
    logging.debug("loading pre-processing for nutation")
    orig_t = s.getaxis("t")
    s.set_units("p_90", "s")
    s.reorder("t", first=True)
    s.chunk("t", ["ph2", "ph1", "t2"], [2, 2, -1])
    s.setaxis("ph2", r_[0.0, 2.0] / 4)
    s.setaxis("ph1", r_[0.0, 2.0] / 4)
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.reorder("t2", first=False)
    s.ft(["ph2", "ph1"])
    if fl is not None:
        fl.next("after phase cycle FT")
        fl.image(s["ph1", 1]["ph2", 0].C.human_units())
    s.ft("t2", shift=True)
    if fl is not None:
        fl.next("freq domain")
        fl.image(s)
    return s


def proc_nutation_amp(s, fl=None):
    logging.debug("loading pre-processing for nutation")
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    if fl is not None:
        fl.next("look for drift")
        fl.image(
            s.C.smoosh(["ph2", "ph1"], "transient")
            .reorder("transient")
            .setaxis("transient", "#")
            .run(abs),
            interpolation="bilinear",
        )
    s.reorder(["ph1", "ph2"])
    s.setaxis("ph2", r_[0:2] / 4).setaxis("ph1", r_[0:4] / 4)
    if "p_90" in s.dimlabels:
        s.set_units("p_90", "s")
    s.ft(["ph1", "ph2"])
    return s


def proc_nutation_chunked(s, fl=None):
    logging.debug("loading pre-processing for nutation")
    s.reorder(["ph1", "ph2"])
    s.set_units("t2", "s")
    s.set_units("p_90", "s")
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft(["ph1", "ph2"])
    s.reorder(["ph1", "ph2", "p_90"])
    if fl is not None:
        fl.next("Raw Data - Time Domain")
        fl.image(s.C.human_units())
    s.ft("t2", shift=False)
    if fl is not None:
        fl.next("Raw Data- Frequency Domain")
        fl.image(s)
    return s


def proc_nutation_v2(s, fl=None):
    logging.debug("loading pre-processing for nutation")
    s.set_units("indirect", "s")
    s.ft(["ph1"])
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.reorder(["ph1"]).set_units("t2", "s")
    if fl is not None:
        fl.next("Raw Data - Time Domain")
        fl.image(s.human_units())
    s.ft("t2", shift=True)
    if fl is not None:
        fl.next("Raw Data- Frequency Domain")
        fl.image(s)
    return s


def proc_nutation_v4(s, fl=None):
    if s.shape["indirect"] > s.shape["nScans"]:
        s.reorder(["ph1", "nScans", "indirect"])
    else:
        s.reorder(["ph1", "indirect", "nScans"])
    s["indirect"] *= 1e-6  # why is it labeled in μs??
    s.rename("indirect", "p_90")
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s.set_units("p_90", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("ph1")
    s.ft("t2", shift=True)
    return s


def proc_var_tau(s, fl=None):
    s.get_prop("SW")
    if "ph1" not in s.dimlabels:
        s.chunk("t", ["ph2", "ph1", "t2"], [2, 4, -1])
        s.setaxis("ph2", r_[0, 2] / 4)
        s.setaxis("ph1", r_[0:4] / 4)
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")  # this should already be set -- why not?
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s *= 2e-6 / 1.11e4  # convert from SpinCore to V (amp)
    s.set_units("V")
    if fl is not None:
        fl.next("raw signal!")
    s.ft("t2", shift=True).ft(["ph1", "ph2"])
    s.reorder(["ph1", "ph2", "tau"])
    if fl is not None:
        fl.plot(abs(s).smoosh(["ph2", "ph1", "tau"], "transients"), alpha=0.2)
        fl.next("raw signal")
        fl.image(s)
    return s


def proc_spincore_echo_v1(s, fl=None):
    "old-fashioned (not properly shaped before storage) echo data"
    s.chunk("t", ["ph2", "ph1", "t2"], [2, 4, -1])
    s.labels({"ph2": r_[0.0, 2.0] / 4, "ph1": r_[0.0, 1.0, 2.0, 3.0] / 4})
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")
    s.reorder("t2", first=False)
    s *= s.shape["nScans"]
    s.squeeze()
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    if "nScans" in s.dimlabels:
        s.setaxis("nScans", "#")
    s.ft("t2", shift=True)
    s.ft(["ph1", "ph2"])
    return s


def proc_spincore_ODNP_v1(s, fl=None):
    logging.debug("loading pre-processing for ODNP")
    prog_power = s.getaxis("power").copy()
    logging.debug(strm("programmed powers", prog_power))
    s.setaxis("power", r_[0 : len(s.getaxis("power"))])
    logging.debug(strm("meter powers", s.get_prop("meter_powers")))
    logging.debug(strm("actual powers", s.getaxis("power")))
    logging.debug(
        strm(
            "ratio of actual to programmed power",
            s.getaxis("power") / prog_power,
        )
    )
    nPoints = s.get_prop("acq_params")["nPoints"]
    SW_kHz = s.get_prop("acq_params")["SW_kHz"]
    nScans = s.get_prop("acq_params")["nScans"]
    nPhaseSteps = s.get_prop("acq_params")["nPhaseSteps"]
    s.chunk("t", ["ph1", "t2"], [4, -1])
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s.labels({"ph1": r_[0.0, 1.0, 2.0, 3.0] / 4})
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft(["ph1"])  # Fourier Transforms coherence channels
    s.reorder(["ph1", "power"])
    s.C.setaxis("power", "#").set_units("power", "scan #")
    if fl is not None:
        fl.next("all data: frequency domain")
        fl.image(s.C.setaxis("power", "#").set_units("power", "scan #"))
    # {{{ since the power axis was saved with settings and not meter powers, fix that here
    p_axis = s.getaxis("power")
    power_axis_dBm = array(s.get_prop("meter_powers"))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = 1e-2 * 10 ** ((power_axis_dBm[:] + 10.0) * 1e-1)
    power_axis_W = r_[0, power_axis_W]
    s.setaxis("power", power_axis_W)
    # }}}
    return s


def proc_spincore_ODNP_v2(s, fl=None):
    logging.debug("loading pre-processing for ODNP")
    prog_power = s.getaxis("power").copy()
    logging.debug(strm("programmed powers", prog_power))
    s.setaxis("power", r_[0 : len(s.getaxis("power"))])
    logging.debug(strm("meter powers", s.get_prop("meter_powers")))
    logging.debug(strm("actual powers", s.getaxis("power")))
    logging.debug(
        strm(
            "ratio of actual to programmed power",
            s.getaxis("power") / prog_power,
        )
    )
    nPoints = s.get_prop("acq_params")["nPoints"]
    SW_kHz = s.get_prop("acq_params")["SW_kHz"]
    nScans = s.get_prop("acq_params")["nScans"]
    nPhaseSteps = s.get_prop("acq_params")["nPhaseSteps"]
    s.chunk("t", ["ph2", "ph1", "t2"], [2, 4, -1])
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")
    s.setaxis("ph2", r_[0.0, 2.0] / 4)
    s.setaxis("ph1", r_[0:4.0] / 4)
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft(["ph1", "ph2"])  # Fourier Transforms coherence channels
    s.C.setaxis("power", "#").set_units("power", "scan #")
    s.reorder(["ph1", "ph2", "power"])
    if fl is not None:
        fl.next("all data: frequency domain")
        fl.image(s.C.setaxis("power", "#").set_units("power", "scan #"))
    # {{{ since the power axis was saved with settings and not meter powers, fix that here
    p_axis = s.getaxis("power")
    power_axis_dBm = array(s.get_prop("meter_powers"))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = 1e-2 * 10 ** ((power_axis_dBm[:] + 10.0) * 1e-1)
    power_axis_W = r_[0, power_axis_W]
    s.setaxis("power", power_axis_W)
    # }}}
    return s


def proc_spincore_ODNP_v3(s, fl=None):
    if "t" in s.dimlabels:
        t.chunk("t", ["ph1", "t2"], [4, -1])
        s.setaxis("ph1", r_[0.0, 1.0, 2.0, 3.0] / 4)
    if "indirect" in s.dimlabels:
        s.rename("indirect", "power")
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s.rename("power", "time")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft(["ph1"])
    if fl is not None:
        fl.next("Raw Data \n Frequency Domain")
        fl.image(s)
        s.ift("t2")
        fl.next("Raw Data \n Time Domain")
        fl.image(s)
        s.ft("t2")
    return s


def proc_spincore_ODNP_v4(s, fl=None):
    if "t" in s.dimlabels:
        t.chunk("t", ["ph2", "ph1", "t2"], [4, 4, -1])
        s.set_units("t2", "s")
    s.rename("power", "time")
    s.setaxis("ph1", r_[0, 1, 2, 3.0] / 4)
    s.setaxis("ph2", r_[0, 1, 2, 3.0] / 4)
    s.set_prop("coherence_pathway", {"ph1": 1, "ph2": -2})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft(["ph1", "ph2"])
    s.reorder(["ph1", "ph2", "time"])
    if fl is not None:
        fl.next("Raw Data \n Frequency Domain")
        fl.image(s)
        s.ift("t2")
        fl.next("Raw Data \n Time Domain")
        fl.image(s)
        s.ft("t2")
    return s


def proc_capture(s):
    logging.debug("loading pre-processing for square wave capture")
    s.set_units("t", "s").name("Amplitude").set_units("V")
    return s


def proc_DOSY_CPMG(s):
    logging.debug("loading pre-processing for DOSY-CPMG")
    # {{{ all of this would be your "preprocessing" and would be tied to the name of your pulse sequence
    l22 = int(s.get_prop("acq")["L"][22])  # b/c the l are integers by definition
    l25 = int(s.get_prop("acq")["L"][25])
    d12 = s.get_prop("acq")["D"][12]
    d11 = s.get_prop("acq")["D"][11]
    p1 = s.get_prop("acq")["P"][1]
    ppg = s.get_prop("pulprog")
    # {{{ these are explanatory -- maybe comment them out?
    m = re.search((".*dwdel1=.*"), ppg, flags=re.IGNORECASE)
    logging.debug(strm(m.groups()))  # show the line that sets dwdel1
    # then look for de and depa
    logging.debug(
        strm(
            [
                (j, s.get_prop("acq")[j])
                for j in s.get_prop("acq").keys()
                if "de" in j.lower()
            ]
        )
    )
    # I actually can't find depa
    # }}}
    m = re.search("\ndefine list<grad_scalar> gl1 = {(.*)}", ppg)
    grad_list = array(
        [float(j.group()) for j in re.finditer("([0-9.]+)", m.groups()[0])]
    )
    m = re.search("([0-9.]+) G/mm", s.get_prop("gradient_calib"))
    grad_list *= float(m.groups()[0]) * 0.1
    dwdel1 = 3.5e-6  # where does this come from? DE is actually larger than this?
    # {{{ find anavpt without hard-setting
    m = re.search('"anavpt=([0-9]+)"', ppg)
    if m is None:
        raise ValueError("I can't find anavpt in the pulse sequence")
    anavpt = int(m.groups()[0])
    # }}}
    dwdel2 = (anavpt * 0.05e-6) / 2
    TD = s.get_prop("acq")["TD2"]
    quadrature_points = TD / 2
    num_points_per_echo = quadrature_points / l25
    acq_time = dwdel2 * num_points_per_echo * 2
    # {{{ so, in principle, later, we can/should do what I did above (w/ eval),
    # but it's getting crazy now, so I stop for now
    tau_extra = 20e-6
    tau_pad = tau_extra - 6e-6
    tau_pad_start = tau_extra - dwdel1 - 6e-6
    tau_pad_end = tau_extra - 6e-6
    tE = (
        dwdel1
        + 5e-6
        + tau_pad_start
        + 1e-6
        + num_points_per_echo * (dwdel2 * 2)
        + tau_pad_end
    )
    # }}}
    s.chunk("indirect", ["indirect", "phcyc"], [l22, -1])
    s.chunk("phcyc", ["ph8", "ph4", "m", "n"], [2, 2, 2, 2])
    s.setaxis("ph8", r_[0.0, 2.0] / 4)
    s.setaxis("ph4", r_[0.0, 2.0] / 4)
    s.setaxis("m", r_[0, 2.0] / 4)
    s.setaxis("n", r_[0, 2.0] / 4)
    s.ft(["ph8", "ph4", "m", "n"])
    s.reorder(["m", "n", "ph4", "ph8", "indirect", "t2"])
    s.setaxis("indirect", grad_list)
    fl.next("abs raw data")
    fl.image(abs(s))
    s.chunk("t2", ["echo", "t2"], [l25, -1])
    s.reorder(["m", "n", "ph4", "ph8", "indirect", "echo", "t2"])
    s.ft("t2", shift=True).ift(
        "t2"
    )  # this is overkill -- need a pyspecdata function that does this w/out the fft
    # }}}
    return s


def proc_ESR(s):
    logging.debug("loading preprocessing for ESR linewidth calculation")
    s -= s["$B_0$", :50].C.mean("$B_0$")
    s_integral = s.C.run_nopop(np.cumsum, "$B_0$")
    x1, x2 = s_integral.getaxis("$B_0$")[r_[5, -5]]
    y1 = s_integral.data[:5].mean()
    y2 = s_integral.data[-5:].mean()
    straight_baseline = (s.fromaxis("$B_0$") - x1) * (y2 - y1) / (x2 - x1)
    s_integral -= straight_baseline
    s_integral /= s_integral.data.mean()
    center_field = (s_integral * s.fromaxis("$B_0$")).mean("$B_0$").item()
    s.setaxis("$B_0$", lambda x: x - center_field)
    s_integral = s.C.run_nopop(np.cumsum, "$B_0$")
    logging.debug(strm(s_integral))
    return s


def proc_field_sweep_v1(s):
    logging.debug(
        "WARNING WARNING, you are using the wrong version of the field sweep code -- should be chunked when data is saved, not on loading!"
    )
    logging.debug("loading preprocessing for fieldsweep")
    s.reorder("t", first=True)
    s.chunk("t", ["ph1", "t2"], [4, -1])
    s.setaxis("ph1", r_[0.0, 1.0, 2.0, 3.0] / 4)
    s.reorder("t2", first=False)
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft("ph1")
    return s


def proc_field_sweep_v2(s):
    s.set_prop("coherence_pathway", {"ph1": 1})
    s.set_units("t2", "s")
    s["t2"] -= s.get_prop("acq_params")["tau_us"] * 1e-6
    s *= s.shape["nScans"]
    s.squeeze()
    s.ft("t2", shift=True)
    s.ft("ph1")
    return s


lookup_table = {
    "chirp": proc_capture,
    "spincore_SE_v1": proc_spincore_SE_v1,
    "spincore_diffph_SE_v1": proc_spincore_diffph_SE_v1,
    "spincore_diffph_SE_v2": proc_spincore_diffph_SE_v2,
    "proc_Hahn_echoph": proc_Hahn_echoph,
    "spincore_IR_v1": proc_spincore_IR,  # for 4 x 2 phase cycle
    "spincore_IR_v2": proc_spincore_IR_v2,  # for 4 x 4 phase cycle data
    "spincore_nutation_v1": proc_nutation,
    "spincore_nutation_v2": proc_nutation_v2,
    "spincore_nutation_amp": proc_nutation_amp,
    "spincore_nutation_v3": proc_nutation_chunked,
    "spincore_nutation_v4": proc_nutation_v4,
    "spincore_ODNP_v1": proc_spincore_ODNP_v1,  # for 4 x 1 phase cycle take meter power
    "spincore_ODNP_v2": proc_spincore_ODNP_v2,  # for 2 x 2 phase cycle take meter powers
    "spincore_ODNP_v3": proc_spincore_ODNP_v3,  # for 4 x 1 phase cycle no meter powers
    "spincore_ODNP_v4": proc_spincore_ODNP_v4,  # for 4 x 4 phase cycle no meter powers
    "spincore_echo_v1": proc_spincore_echo_v1,
    "spincore_var_tau_v1": proc_var_tau,
    "square_wave_capture_v1": proc_capture,
    "DOSY_CPMG_v1": proc_DOSY_CPMG,
    "ESR_linewidth": proc_ESR,
    "field_sweep_v1": proc_field_sweep_v1,
    "field_sweep_v2": proc_field_sweep_v2,
}
