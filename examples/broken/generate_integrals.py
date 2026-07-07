"""
Convert 2D to Integrals with Errors
===================================

Take a 2D dataset and convert it to a table of integrals with errors, utilizing
all the bells and whistles (frequency and time selection, alignment, etc.)

Demonstrate on a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise, and with fluctuating field
(normally distributed field variation).
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from numpy.random import seed
import sympy as s
from collections import OrderedDict

init_logging(level="debug")

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 3

fl = fl_mod()
t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
echo_time = 10e-3
with figlist_var() as fl:
    for (
        expression,
        variable_defs,
        signal_pathway,
        indirect,
        clock_correction,
        label,
        f_range,
    ) in [
        (
            (
                23
                * (1 - 2 * s.exp(-vd / 0.2))
                * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi)
            ),
            [
                ("vd", nddata(r_[0:1:40j], "vd")),
                ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
                ("ph2", nddata(r_[0, 2] / 4.0, "ph2")),
                ("t2", nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 0, "ph2": 1},
            "vd",
            True,
            "IR",
            (-400, 400),
        ),
        (
            (
                23
                * (1 - (3.2 * power / (0.25 + power)) * 150e-6 * 659.33)
                * s.exp(+1j * 2 * s.pi * 100 * t2 - abs(t2) * 50 * s.pi)
            ),
            [
                ("power", nddata(r_[0:4:25j], "power")),
                ("ph1", nddata(r_[0:4] / 4.0, "ph1")),
                ("t2", nddata(r_[0:0.2:256j] - echo_time, "t2")),
            ],
            {"ph1": 1},
            "power",
            False,
            "Enhancement",
            (-200, 600),
        ),
    ]:
        fl.basename = "(%s)" % label
        data = fake_data(
            expression, OrderedDict(variable_defs), signal_pathway
        )
        data.ft("t2")
        # {{{ make data unitary again
        data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop("t2", "dt")
        # }}}
        zero_pathway = {j: 0 for j in signal_pathway}
        excluded_pathways = [signal_pathway]
        if zero_pathway != signal_pathway:
            excluded_pathways.append(zero_pathway)
        if clock_correction and "nScans" in data.dimlabels:
            data = clock_correct(data, indirect=indirect, fl=fl)
        data = data["t2":f_range]
        # TODO ☐: in the following, you need
        # to investigate the git history to
        # explain where this was removed, and
        # summarize the result of any added
        # comments or git commit msgs.  This
        # is also relevant to why we need to
        # call clock correct above
        # The old generate_integrals helper is gone; this example now
        # calls the lower-level routine that performs its key final step:
        # find frequency bounds, integrate the selected pathway, and
        # propagate the error.
        data_int = frequency_domain_integral(
            data,
            signal_pathway=signal_pathway,
            excluded_pathways=excluded_pathways,
            indirect=indirect,
            fl=fl,
        )
        fl.next(f"{label} integrals")
        fl.plot(data_int, ".", capsize=6)
