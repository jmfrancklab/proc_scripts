"""QESR of same sample
===================

Here, we sit and acquire the same sample with different acquisition
parameters, and make sure we get the same QESR.

We do this b/c the data is not stored by XEPR in a way that we would
expect.

Note that some of these datasets are not great b/c we didn't capture
enough of the baseline to either side (and we see here why this is a
problem), and so the QESR doesn't do great.

However, both with the single integral (absorption) plot and more
precisely with the ESR alignment we conduct at the end,
we see that we have properly scaled the spectra.
Note this doesn't mean they are equal -- we are able to
see subtle (or not so subtle) differences in the
scaling because of overmodulation -- if the spectra are
not overmodulated, they match in scaling to within a
few percent.

See the `QESR.py` example for information about setting your pyspecdata
config so that this works correctly!
"""

from matplotlib.pyplot import title, rcParams
from pyspecdata import find_file, figlist_var
from pyspecProcScripts import QESR, align_esr
from collections import OrderedDict

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 9

fieldaxis = "$B_0$"
plot_rescaled = False
water_nddata = find_file("230511_water.DSC", exp_type="francklab_esr/romana")[
    "harmonic", 0
]
d = OrderedDict()
for file_searchstring, thislabel, exp_type in [
    ("250321_OHT_control", "control", "francklab_esr/Warren"),
    ("250321_OHT_gaindown", "lower gain", "francklab_esr/Warren"),
    ("250321_OHT_moddown", "lower mod", "francklab_esr/Warren"),
    (
        "250321_OHT_followmodrule",
        "lowest (good) mod",
        "francklab_esr/Warren",
    ),
    ("250321_OHT_morepoint", "more points", "francklab_esr/Warren"),
    ("250321_OHT_lowpower", "low power", "francklab_esr/Warren"),
]:
    d[thislabel] = find_file(file_searchstring, exp_type=exp_type).chunk_auto(
        "harmonic"
    )["harmonic", 0]["phase", 0]
with figlist_var() as fl:
    fl.next("Raw and unscaled")
    for k, v in d.items():
        fl.plot(v, label=k, alpha=0.7)
    for k, v in d.items():
        c = QESR(
            v,
            label=k,
            fl=fl,
        )
        print(f"{c:~P}")
    d = align_esr(d, fl=fl)
    title(
        "This plot shows we know how to scale the amplitude\n"
        "correctly for a range of different experimental parameters"
    )
