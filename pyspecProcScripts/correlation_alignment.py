import pyspecdata as psd
import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import logging


@FuncFormatter
def to_percent(y, position):
    s = "%.2g" % (100 * y)
    if plt.rcParams["text.usetex"] is True:
        return s + r"$\%$"
    else:
        return s + "%"


def correl_align(
    s_jk,
    align_phases=False,
    tol=1e-4,
    indirect_dim="indirect",
    fig_title="correlation alignment",
    signal_pathway={"ph1": 0, "ph2": 1},
    shift_bounds=False,
    avg_dim=None,
    max_shift=100.0,
    sigma=20.0,
    direct="t2",
    fl=None,
):
    """
    Align transients collected with chunked phase cycling dimensions along an
    indirect dimension based on maximizing the correlation across all the
    transients and repeat alignment until the calculated signal energy remains
    constant to within a given tolerance level.

    Parameters
    ==========
    s_jk:  psd.nddata
        A psd.nddata object which contains phase cycle dimensions and an
        indirect dimension. This data corresponds to the s_jk signal in the
        correlation equations in the DCCT paper
    align_phases:   boolean
    tol:            float
                    Sets the tolerance limit for the alignment procedure.
    indirect_dim:   str
                    Name of the indirect dimension along which you seek to
                    align
                    the transients.
    fig_title:      str
                    Title for the figures generated.
    signal_pathway: dict
                    Dictionary containing the signal pathway.
    shift_bounds:   boolean
                    Keeps f_shift to be within a specified
                    limit (upper and lower bounds given by max_shift)
                    which should be around the location of the expected
                    signal.
    avg_dim:        str
                    Dimension along which the data is being averaged.
    max_shift:      float
                    Specifies the upper and lower bounds to the range over
                    which f_shift will be taken from the correlation function.
                    Shift_bounds must be True.
    sigma:          int
                    Sigma value for the Gaussian mask. Related to the line
                    width of the given data.
    fl:             boolean
                    fl=fl to show the plots and figures produced by this
                    function otherwise, fl=None.

    Returns
    =======
    f_shift:    np.array
                The optimized frequency shifts for each transient which will
                maximize their correlation amongst each other, thereby aligning
                them.
    sigma:      float
                The width of the Gaussian function used to frequency filter
                the data in the calculation of the correlation function.
    """

    logging.debug(psd.strm("Applying the correlation routine"))
    if avg_dim:
        phcycdims = [j for j in s_jk.dimlabels if j.startswith("ph")]
        indirect = set(s_jk.dimlabels) - set(phcycdims) - set([direct])
        indirect = [j for j in s_jk.dimlabels if j in indirect]
        avg_dim_len = len(s_jk.getaxis(avg_dim))
        s_jk.smoosh(indirect)
    for unk_key in signal_pathway.keys():
        assert not s_jk.get_ft_prop(unk_key), (
            str(unk_key) + " must not be in the coherence domain"
        )
    # {{{ Determine names and values of CT pathways and indirect    
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    ph_len_dict = {j: psd.ndshape(s_jk)[j] for j in signal_pathway.keys()}
    N = psd.ndshape(s_jk)[indirect_dim]
    # }}}
    # energy of unaligned signal
    sig_energy = (abs(s_jk) ** 2).data.sum().item() / N
    # {{{ Diagnostic plot for unaligned
    if fl:
        fl.push_marker()
        fig_forlist, ax_list = plt.subplots(1, 5, figsize=(25, 10))
        fl.next("Correlation Diagnostics")
        fig_forlist.suptitle(
            " ".join(
                ["Correlation Diagnostic"]
                + [j for j in [fl.basename] if j is not None]
            )
        )
        s_jk.reorder([direct], first=False)
        fl.image(
            s_jk.C.setaxis(indirect_dim, "#").set_units(
                indirect_dim, "scan #"
            ),
            ax=ax_list[0],
            human_units=False,
        )
        ax_list[0].set_title("before correlation\nsig. energy=%g" % sig_energy)
    # }}}    
    energy_diff = 1.0
    i = 0
    energy_vals = []
    this_E = (abs(s_jk.C.sum(indirect_dim)) ** 2).data.sum().item() / N**2
    energy_vals.append(this_E / sig_energy)
    last_E = None
    # {{{ find center frequency to see where to center the mask
    for_nu_center = s_jk.C
    for_nu_center.ft(list(signal_pathway))
    for x in range(len(signal_keys)):
        for_nu_center = for_nu_center[signal_keys[x], signal_values[x]]
    nu_center = for_nu_center.mean(indirect_dim).C.argmax(direct)
    # }}}
    logging.debug(psd.strm("Center frequency", nu_center))
    for my_iter in range(100):
        i += 1
        print("HEY",i)
        logging.debug(psd.strm("*** *** ***"))
        logging.debug(psd.strm("CORRELATION ALIGNMENT ITERATION NO. ", i))
        logging.debug(psd.strm("*** *** ***"))
        s_jk.ift(direct)
        # s_mn will become the 3 dimensional signal with the additional
        # Delta ph# dimension. This corresponds to the signal s_mn in the
        # DCCT paper. But because we need the original unaligned later for
        # diagnonstics we make a copy
        s_mn = s_jk.C
        if align_phases:
            ph0 = s_jk.C.sum(direct)
            ph0 /= abs(ph0)
            s_mn /= ph0
        # {{{ Apply mask around center of signal in nu domain
        s_mn.ft(direct)
        this_mask = np.exp(
            -((s_mn.fromaxis(direct) - nu_center) ** 2) / (2 * sigma**2)
        )
        s_mn *= this_mask
        # }}}
        s_mn.ift(direct)
        # {{{ Make extra dimension for s_mn for the 
        for ph_name, ph_len in ph_len_dict.items():
            ph = np.ones(ph_len)
            s_mn *= psd.nddata(ph, "Delta" + ph_name.capitalize())
            s_mn.setaxis("Delta" + ph_name.capitalize(), "#")
        correl = s_mn * 0
        # {{{ phase correction term in eq 29 in DCCT paper
        for ph_name, ph_len in ph_len_dict.items():
            for ph_index in range(ph_len):
                s_mn["Delta%s" % ph_name.capitalize(), ph_index] = s_mn[
                    "Delta%s" % ph_name.capitalize(), ph_index
                ].run(lambda x, axis=None: np.roll(x, ph_index, axis=axis), ph_name)
        # }}}        
        correl = s_mn.mean(indirect_dim).run(np.conj) * s_jk
        correl.reorder([indirect_dim, direct], first=False)
        if my_iter == 0:
            logging.debug(psd.strm("holder"))
            if fl:
                correl.reorder([direct], first=False)
                fl.image(
                    correl.C.setaxis(indirect_dim, "#").set_units(
                        indirect_dim, "scan #"
                    ),
                    ax=ax_list[1],
                )
                ax_list[1].set_title("correlation function (t), \nafter apod")
        # we want Δnu centered at 0 and don't care about the previous window selection
        correl.ft_new_startpoint(direct, "freq")
        correl.ft_new_startpoint(direct, "time")
        # Part of convolution is multiplying in t domain and then zero-filling
        # before re-FT
        correl.ft(direct, shift=True, pad=2**14)
        # {{{ Isolate the correlation function for the signal pathway
        for sig_pathname, sig_pathval in signal_pathway.items():
            correl.ft(["Delta%s" % sig_pathname.capitalize()])
            correl = (
                correl["Delta" + sig_pathname.capitalize(), sig_pathval]
                + correl["Delta" + sig_pathname.capitalize(), 0]
            )
        # }}}
        if my_iter == 0:
            logging.debug(psd.strm("holder"))
            if fl:
                correl.reorder([direct], first=False)
                fl.image(
                    correl.C.setaxis(indirect_dim, "#").set_units(
                        indirect_dim, "scan #"
                    ),
                    ax=ax_list[2],
                    human_units=False,
                )
                ax_list[2].set_title("correlation function (v), \nafter apod")
        # Find optimal f shift based on max of correlation function        
        if shift_bounds:
            f_shift = (
                correl[direct:(-max_shift, max_shift)]
                .run(np.real)
                .argmax(direct)
            )
        else:
            f_shift = correl.run(np.real).argmax(direct)
        s_copy = s_jk
        s_copy *= np.exp(-1j * 2 * np.pi * f_shift * s_copy.fromaxis(direct))
        s_jk.ft(direct)
        s_copy.ft(direct)
        if my_iter == 0:
            logging.debug(psd.strm("holder"))
            if fl:
                s_copy.reorder([direct], first=False)
                fl.image(
                    s_copy.C.setaxis(indirect_dim, "#").set_units(
                        indirect_dim, "scan #"
                    ),
                    ax=ax_list[3],
                    human_units=False,
                )
                ax_list[3].set_title("after correlation\nbefore ph0 restore")
        logging.debug(
            psd.strm(
                "signal energy per transient (recalc to check that it stays"
                " the same):",
                (abs(s_copy**2).data.sum().item() / N),
            )
        )

        this_E = (
            abs(s_copy.C.sum(indirect_dim)) ** 2
        ).data.sum().item() / N**2
        energy_vals.append(this_E / sig_energy)
        logging.debug(
            psd.strm("averaged signal energy (per transient):", this_E)
        )
        if last_E is not None:
            energy_diff = (this_E - last_E) / sig_energy
            logging.debug(psd.strm(energy_diff))
            if abs(energy_diff) < tol and my_iter > 4:
                break
        last_E = this_E
        print(f_shift)
    if fl is not None:
        fl.next("correlation convergence")
        fl.plot(np.array(energy_vals), "x")
        plt.gca().yaxis.set_major_formatter(to_percent)
    if fl is not None:
        fl.image(
            s_copy.C.setaxis(indirect_dim, "#").set_units(
                indirect_dim, "scan #"
            ),
            ax=ax_list[4],
        )
        ax_list[4].set_title(
            "after correlation\nph0 restored \nsig. energy=%g" % sig_energy
        )
        fl.pop_marker()
    return f_shift, sigma, this_mask
