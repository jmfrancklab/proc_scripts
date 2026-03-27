from pyspecdata import find_file, DCCT

# from pyspecProcScripts import lookup_table
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyspecProcScripts as prscr
import pyspecdata as psd

plt.rcParams["figure.figsize"] = (10, 10)

apply_apod = False
center_echo = False
freq_filtering = True
exp_list = [
    ("260323", (-7500.0, -2370.0), 3),
    ("260324", (-3330.0, 250.0), 4),
    ("260324", (-10000.0, -6000.0), 6),
]
this_lim = 0
for thisdate, frq_lims, whichexp in exp_list:
    this_node = f"shim_y_{whichexp}"
    d = find_file(
        rf"{thisdate}.*{whichexp}\.h5",
        exp_type="B27/Echoes",
        expno=this_node,
        lookup=prscr.lookup_table,
    )
    beta = d.get_prop("acq_params")["beta_90_s_sqrtW"]
    freq = d.get_prop("acq_params")["carrierFreq_MHz"]
    nscans = d.get_prop("acq_params")["nScans"]
    if nscans > 1:
        d = d.mean("nScans")
    if freq_filtering:
        d = d["t2":frq_lims]
    d.ift("t2")
    # {{{ first figure covers data selection
    fig = plt.figure(f"{this_node} data freq_filtering")
    gs = mpl.gridspec.GridSpec(1, 2, left=0.1, right=0.95, top=0.9, bottom=0.2)
    DCCT(
        d,
        fig=fig,
        title=(
            f"{this_node}\n β={beta} time domain \n f = {freq}",
            f" nscans = {nscans}",
        ),
        bbox=gs[0, 0],
    )
    d = d["t2", 1:]  # remove first point, which has ringdown
    if apply_apod:
        d *= np.exp(
            -abs(d.fromaxis("t2")) / 0.5e-3
        )  # apply a filter to clean up signal
    d.ft("t2")
    DCCT(
        d,
        fig=fig,
        title=(
            f"{this_node}\n β={beta} frequency domain \n f = {freq},"
            f" nscans = {nscans}"
        ),
        interpolation="nearest",
        bbox=gs[0, 1],
    )
    plt.tight_layout()
    # }}}
    # {{{ carry through to energy
    fig = plt.figure(f"{this_node} determine energy")
    gs = mpl.gridspec.GridSpec(
        1, 3, left=0.1, right=0.95, top=0.9, bottom=0.15
    )
    DCCT(
        d,
        fig=fig,
        title=(
            f"{this_node}\n β={beta} frequency domain \n"
            f"f = {freq} nscans = {nscans}"
        ),
        interpolation="nearest",
        bbox=gs[0, 0],
    )
    d = abs(prscr.select_pathway(d, d.get_prop("coherence_pathway")))
    d.pcolor(ax=fig.add_subplot(gs[0, 1]))
    d.run(lambda x: abs(x) ** 2).integrate("t2")
    plt.tight_layout()
    psd.plot(d, "o", ax=fig.add_subplot(gs[0, 2]))
    plt.title(this_node)
    plt.tight_layout()
    this_lim += 1
    # }}}
plt.show()
