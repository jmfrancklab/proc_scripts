from pyspecdata import find_file, image, DCCT, gammabar_H

# from pyspecProcScripts import lookup_table
import matplotlib.pyplot as plt
from numpy import r_, pi
import matplotlib as mpl
import numpy as np
import pyspecProcScripts as prscr
from pyspecProcScripts.load_data import proc_spincore_generalproc_v1
import pyspecdata as psd

plt.rcParams["figure.figsize"] = (10, 10)

apply_apod = False
center_echo = False
slicing = True
lims = (-7.1e3, -1.5e3)

for whichexp in [2, 3]:
    this_expno = f"shim_y_{whichexp}"
    d = find_file(
        rf"260323.*{whichexp}\.",
        exp_type="B27/Echoes",
        expno=this_expno,  # when I leave out expno, it gives me an error that lists all the nodes
        lookup=prscr.lookup_table,
    )
    beta = d.get_prop("acq_params")["beta_90_s_sqrtW"]
    freq = d.get_prop("acq_params")["carrierFreq_MHz"]
    ind_pts = d.get_prop("acq_params")["indirect_pts"]
    d.ift("t2")
    # {{{ first figure covers data selection
    fig = plt.figure(f"{this_expno} data slicing")
    gs = mpl.gridspec.GridSpec(1, 2, left=0.1, right=0.95, top=0.9, bottom=0.2)
    DCCT(
        d,
        fig=fig,
        title=(
            f"{this_expno}\n β={beta} time domain \n f = {freq} ind_pts = {ind_pts}"
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
            f"{this_expno}\n β={beta} frequency domain \n f = {freq} ind_pts ="
            f" {ind_pts}"
        ),
        interpolation="nearest",
        bbox=gs[0, 1],
    )
    plt.tight_layout()
    # }}}
    # {{{ carry through to energy
    fig = plt.figure(f"{this_expno} determine energy")
    gs = mpl.gridspec.GridSpec(1, 3, left=0.1, right=0.95, top=0.9, bottom=0.15)
    if slicing:
        d = d["t2":lims]
    DCCT(
        d,
        fig=fig,
        title=(
            f"{this_expno}\n β={beta} frequency domain \n f = {freq} ind_pts ="
            f" {ind_pts}"
        ),
        interpolation="nearest",
        bbox=gs[0, 0],
    )
    d = abs(prscr.select_pathway(d, d.get_prop("coherence_pathway")))
    d.pcolor(ax=fig.add_subplot(gs[0, 1]))
    d.run(lambda x: abs(x) ** 2).integrate("t2")
    plt.tight_layout()
    psd.plot(d, "o", ax=fig.add_subplot(gs[0, 2]))
    plt.title(this_expno)
    plt.tight_layout()
    # }}}
plt.show()
