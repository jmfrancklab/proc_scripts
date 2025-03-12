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
    s_orig,
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

    This implements Eq 29 in the DCCT paper (Beaton et al JCP 2022),
    remembering that the FT of the cross-correlation is the product
    of the two functions, with the *first* on complex conjugate.

    Parameters
    ==========
    s_orig:  psd.nddata
        A psd.nddata object which contains phase cycle dimensions and an
        indirect dimension.
        This data is not modified.
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
    # TODO ☐: I think we want to remove avg_dim, but let's wait and see.
    # TODO ☐: I think that we want to just allow the following to
    #         determine the indirect dimensions.
    # TODO ☐: I think we want to rename and throw errors for the kwargs
    #         "indirect_dim" and "avg_dim", and rather use kwargs that
    #         clearly specify which indirect dimensions are *used for
    #         alignment* and which are not.
    #         For example, in an inversion recovery or E(p) we typically
    #         want to align along the indirect dimension, but only AFTER
    #         we have made it safe to do so by flipping the signs.
    #         As another example, we almost *always* want to align along
    #         nScans.
    if avg_dim:
        phcycdims = [j for j in s_orig.dimlabels if j.startswith("ph")]
        indirect = set(s_orig.dimlabels) - set(phcycdims) - set([direct])
        indirect = [j for j in s_orig.dimlabels if j in indirect]
        s_jk = s_orig.C.smoosh(indirect)  # this version ends up with
        #                                   three dimensions
        #                                   (j=align_dim, k=phcyc, and
        #                                   direct nu) and is NOT conj
    else:
        s_jk = s_orig.C  # even if there isn't an indirect to smoosh we will
        #                 later be applying modifications to s_jk that we don't
        #                 want applied to s_orig
    for phnames in signal_pathway.keys():
        assert not s_orig.get_ft_prop(phnames), (
            str(phnames) + " must not be in the coherence domain"
        )
    assert s_orig.get_ft_prop(direct), "direct dimension must be in the frequency domain"
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    ph_len = {j: psd.ndshape(s_orig)[j] for j in signal_pathway.keys()}
    N = psd.ndshape(s_orig)[indirect_dim]
    sig_energy = (abs(s_orig) ** 2).data.sum().item() / N
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
        # TODO ☐: why isn't this a DCCT?
        fl.image(
            s_jk,
            ax=ax_list[0],
            human_units=False,
        )
        ax_list[0].set_title("before correlation\nsig. energy=%g" % sig_energy)
    energy_diff = 1.0
    i = 0
    energy_vals = []
    this_E = (abs(s_orig.C.sum(indirect_dim)) ** 2).data.sum().item() / N**2
    energy_vals.append(this_E / sig_energy)
    last_E = None
    # TODO ☐: the mask needs to be separated into a different function
    # {{{ find center frequency to see where to center the mask
    # TODO ☐: this copy is undesirable, but not dealing with it, since
    #         we need to separate the mask anyways
    for_nu_center = s_jk.C
    for_nu_center.ft(list(signal_pathway))
    for x in range(len(signal_keys)):
        for_nu_center = for_nu_center[signal_keys[x], signal_values[x]]
    nu_center = for_nu_center.mean(indirect_dim).C.argmax(direct)
    logging.debug(psd.strm("Center frequency", nu_center))
    # }}}
    # {{{ Apply mask around center of signal in frequency domain.
    #     This will become the expression in the left brackets, where
    #     s_mn is at this stage equal to s_jk.
    this_mask = np.exp(
        -((s_jk.fromaxis(direct) - nu_center) ** 2) / (2 * sigma**2)
    )
    s_leftbracket = this_mask * s_jk
    # }}}
    s_leftbracket.ift(direct)
    s_jk.ift(direct)
    for my_iter in range(100):
        i += 1
        logging.debug(psd.strm("*** *** ***"))
        logging.debug(psd.strm("CORRELATION ALIGNMENT ITERATION NO. ", i))
        logging.debug(psd.strm("*** *** ***"))
        # {{{ optionally make the phases of all the different phase
        #     cycle steps the same
        if align_phases:
            ph0 = s_jk.C.sum(direct)
            ph0 /= abs(ph0)
            s_jk /= ph0
        # }}}
        # {{{ Make extra dimension (Δφ) for s_leftbracket
        for phname, phlen in ph_len.items():
            ph = np.ones(phlen)
            s_leftbracket *= psd.nddata(ph, "Delta" + phname.capitalize())
            s_leftbracket.setaxis("Delta" + phname.capitalize(), "#")
        # }}}
        correl = s_leftbracket * 0
        # {{{ phase correction term in eq 29 in DCCT paper
        for phname, phlen in ph_len.items():
            for ph_index in range(phlen):
                s_leftbracket["Delta%s" % phname.capitalize(), ph_index] = s_leftbracket[
                    "Delta%s" % phname.capitalize(), ph_index
                ].run(
                    lambda x, axis=None: np.roll(x, ph_index, axis=axis),
                    phname,
                )
        # }}}
        correl = s_leftbracket.mean(indirect_dim).run(np.conj) * s_jk
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
        # we want Δnu centered at 0 and don't care about the previous window
        # selection
        correl.ft_new_startpoint(direct, "freq")
        correl.ft_new_startpoint(direct, "time")
        # Part of convolution is multiplying in t domain and then zero-filling
        # before re-FT
        correl.ft(direct, shift=True, pad=2**14)
        # {{{ Apply mask and sum along Δp_l
        for ph_name, ph_val in signal_pathway.items():
            # FT transforms from Δφ to Δp_l
            correl.ft(["Delta%s" % ph_name.capitalize()])
            correl = (
                correl["Delta" + ph_name.capitalize(), ph_val]
                + correl["Delta" + ph_name.capitalize(), 0]
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
        # Calculate energy with shift applied
        s_aligned = s_jk * np.exp(
            -1j * 2 * np.pi * f_shift * s_jk.fromaxis(direct)
        )
        s_aligned.ft(direct)
        if my_iter == 0:
            logging.debug(psd.strm("holder"))
            if fl:
                fl.image(
                    s_orig.C.setaxis(indirect_dim, "#").set_units(
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
                (abs(s_aligned**2).data.sum().item() / N),
            )
        )
        # {{{ Calculate energy difference from last shift to see if
        #     there is any further gain to keep reiterating
        this_E = (
            abs(s_aligned.C.sum(indirect_dim)) ** 2
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
        # }}}
        last_E = this_E
    if fl is not None:
        fl.next("correlation convergence")
        fl.plot(np.array(energy_vals), "x")
        plt.gca().yaxis.set_major_formatter(to_percent)
    if fl is not None:
        fl.image(
            s_aligned.C.setaxis(indirect_dim, "#").set_units(
                indirect_dim, "scan #"
            ),
            ax=ax_list[4],
        )
        ax_list[4].set_title(
            "after correlation\nph0 restored \nsig. energy=%g" % sig_energy
        )
        fl.pop_marker()
    return f_shift, sigma, this_mask
