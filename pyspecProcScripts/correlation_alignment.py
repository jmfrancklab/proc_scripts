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
    tol=1e-4,
    indirect_dim=None,
    repeat_dims=[],
    non_repeat_dims=[],
    fig_title="correlation alignment",
    signal_pathway={"ph1": 0, "ph2": 1},
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
    tol:            float
                    Sets the tolerance limit for the alignment procedure.
    repeat_dims:   list (default [])
                    List of the dimensions along which the signal is
                    essentially repeated and therefore can be aligned.
                    If there is an nScans or repeats dimension this must
                    be included! By default, this list is empty. Note: prior to
                    calling this function it is a common strategy to flip the
                    sign so all transients are the same sign.
    nonrepeat_dims: list (default [])
                    These are indirect dimensions that we are not OK to align
                    (e.g. the indirect dimension of a 2D COSY experiment).
    fig_title:      str
                    Title for the figures generated.
    signal_pathway: dict
                    Dictionary containing the signal pathway.
    max_shift:      float
                    Specifies the upper and lower bounds to the range over
                    which f_shift will be taken from the correlation function.
                    Shift_bounds must be True.
                    If it's set to None, then no bounds are applied.
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
    logging.debug("Applying the correlation routine")
    assert indirect_dim is None, (
        "We updated the correlation function to no longer take indirect_dim as"
        " a kwarg! Now indirect_dim roughly corresponds to repeat_dims"
    )
    assert avg_dim is None, (
        "We updated the correlation function to no longer take avg_dim as a"
        " kwarg!"
    )
    assert (
        type(repeat_dims) == list
    ), "the repeat_dims kwarg needs to be a list of strings"
    assert (
        len(repeat_dims) > 0
    ), "You must tell me which dimension contains the repeats!"
    phcycdims = [j for j in s_orig.dimlabels if j.startswith("ph")]
    repeats = set(s_orig.dimlabels) - set(phcycdims) - set([direct])
    assert len(repeats - set(repeat_dims) - set(non_repeat_dims)) == 0
    repeat_dims = [j for j in s_orig.dimlabels if j in repeats]
    if len(repeat_dims) > 1:
        s_jk = s_orig.C.smoosh(repeats, "repeats")  # this version ends up with
        #                                   three dimensions
        #                                   (j=align_dim, k=phcyc, and
        #                                   direct nu) and is NOT conj
    else:
        s_jk = s_orig.C.rename(
            repeat_dims, "repeats"
        )  # even if there isn't an indirect to smoosh we will
        #                 later be applying modifications to s_jk that we don't
        #                 want applied to s_orig
    for phnames in signal_pathway.keys():
        assert not s_orig.get_ft_prop(phnames), (
            str(phnames) + " must not be in the coherence domain"
        )
    assert s_orig.get_ft_prop(
        direct
    ), "direct dimension must be in the frequency domain"
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    ph_len = {j: psd.ndshape(s_orig)[j] for j in signal_pathway.keys()}
    N = psd.ndshape(s_jk)["repeats"]
    # TODO ☐: as noted below, this doesn't include the mask!
    sig_energy = (abs(s_jk) ** 2).data.sum().item() / N
    if fl:
        fl.push_marker()
        fig_forlist, ax_list = plt.subplots(1, 4, figsize=(25, 10))
        fl.next("Correlation Diagnostics")
        fig_forlist.suptitle(
            " ".join(
                ["Correlation Diagnostic"]
                + [j for j in [fl.basename] if j is not None]
            )
        )
        s_jk.reorder([direct], first=False)
        # TODO ☐: why isn't this a DCCT? (and other calls to image in
        #         this function)
        fl.image(
            s_jk,
            ax=ax_list[0],
            human_units=False,
        )
        ax_list[0].set_title("before correlation\nsig. energy=%g" % sig_energy)
    energy_diff = 1.0
    i = 0
    energy_vals = []
    # below we divide by an extra N because if the signals along the repeats
    # are the same, then the energy of the resulting sum should increase by N
    # (vs taking the square and summing which is what we do for calculating the
    # sig_energy above)
    this_E = (abs(s_jk.C.sum("repeats")) ** 2).data.sum().item() / N**2
    energy_vals.append(this_E / sig_energy)
    last_E = None
    # TODO ☐: the mask needs to be separated into
    #         two functions:
    #         One (f_mask) takes in s_jk and returns a copy
    #         with the *frequency* mask applied.
    #         Another (Delta_p_mask) takes in
    #         s_leftbracket, and returns the
    #         result of applying the mask along
    #         Δp_l
    # {{{ find center frequency to see where to center the mask
    # TODO ☐: this copy is undesirable, but not dealing with it, since
    #         we need to separate the mask
    #         anyways.  Likely, in the final
    #         version, when we supply the mask
    #         function, this will be determined
    #         from the same code that applies the
    #         frequency bounds.
    for_nu_center = s_jk.C
    for_nu_center.ft(list(signal_pathway))
    for x in range(len(signal_keys)):
        for_nu_center = for_nu_center[signal_keys[x], signal_values[x]]
    nu_center = for_nu_center.mean("repeats").C.argmax(direct)
    logging.debug(psd.strm("Center frequency", nu_center))
    # }}}
    s_jk.ift(direct)
    f_shift = 0
    for my_iter in range(100):
        # Note that both s_jk and s_leftbracket
        # change every iteration, because the
        # *data* is updated with every iteration
        i += 1
        logging.debug(psd.strm("*** *** ***"))
        logging.debug(psd.strm("CORRELATION ALIGNMENT ITERATION NO. ", i))
        logging.debug(psd.strm("*** *** ***"))
        # {{{ construct the expression in the left square brackets of
        #     eq. 29.
        # {{{ Apply mask around center of signal
        #     in frequency domain.
        #     At this stage, s_mn is equal to
        #     s_jk.
        #
        #     Note that is seems expensive, but
        #     the masked data does genuinely
        #     change with every iteration, because
        #     the signal frequency is moving
        #     relative to the mask.
        s_jk.ft(direct)
        this_mask = np.exp(
            -((s_jk.fromaxis(direct) - nu_center) ** 2) / (2 * sigma**2)
        )
        s_leftbracket = this_mask * s_jk
        s_jk.ift(direct)
        s_leftbracket.ift(direct)
        # }}}
        # {{{ Make extra dimension (Δφ_n) for s_leftbracket:
        #     start by simply replicating the data along the new
        #     dimension.
        for phname, phlen in ph_len.items():
            ph_ones = np.ones(phlen)
            s_leftbracket *= psd.nddata(ph_ones, "Delta" + phname.capitalize())
            s_leftbracket.setaxis("Delta" + phname.capitalize(), "#")
        # }}}
        # {{{ Currently, we have a Δφ_n dimension, but we've just
        #     copied/replicated the data, so the phases are like this:
        #
        #     Δφ_n →
        # φ_k 0 0 0 0
        #  ↓  1 1 1 1
        #     2 2 2 2
        #     3 3 3 3
        #     For this to actually correspond to a phase difference
        #     (Δφ_n), we need to roll the data for each "column" of
        #     Δφ_n by an increasing amount, so it looks like this
        #
        #     Δφ_n →
        # φ_k 0 4 3 2
        #  ↓  1 0 4 3
        #     2 1 0 4
        #     3 2 1 0
        for phname, phlen in ph_len.items():
            for ph_index in range(phlen):
                s_leftbracket[
                    "Delta%s" % phname.capitalize(), ph_index
                ] = s_leftbracket[
                    "Delta%s" % phname.capitalize(), ph_index
                ].run(
                    lambda x, axis=None: np.roll(x, ph_index, axis=axis),
                    phname,
                )
        # }}}
        # }}}
        # the sum over m in eq. 29 only applies to the left bracket,
        # so we just do it here
        correl = s_leftbracket.mean("repeats").run(np.conj) * s_jk
        correl.reorder(["repeats", direct], first=False)
        if my_iter == 0:
            logging.debug(psd.strm("holder"))
            if fl:
                correl.reorder([direct], first=False)
                fl.image(
                    correl,
                    ax=ax_list[1],
                )
                ax_list[1].set_title(
                    "correlation function (t)\n(includes ν mask)"
                )
        # {{{ FT the correlation function so that we can determine the
        #     relative shift needed to line each transient up with the
        #     average:
        #     - we need to clear the frequency startpoint and the do a
        #       shifted FT in order to get frequencies that are just
        #       centered at 0 (no shift)
        #     - we need to pad (zero fill time domain) so that we get a
        #       very precise shift value, and are not limited by our
        #       frequency resolution
        correl.ft_new_startpoint(direct, "freq")
        correl.ft_new_startpoint(direct, "time")
        correl.ft(direct, shift=True, pad=2**14)
        # }}}
        # TODO ☐: only the left square bracket term depends on Δp_l,
        #         so the following should be done on the left square
        #         bracket term before multiplication
        # TODO ☐: then, the following sum is
        #         actually part of the masking
        #         procedure, so it should be moved
        #         into the Delta_p_mask that is
        #         supplied here.
        #         The point of this is that right
        #         now, we are ALWAYS maximizing
        #         signal in both the coherence
        #         pathway of interest and the 0,0
        #         pathway → we would like to
        #         control that.
        # {{{ this applies the Fourier transform from Δφ to Δp_l
        #     that is found inside the left square bracket of eq. 29.
        #     Then, it performs a sum that is equivalent to applying a
        #     mask along the Δp_l dimension and then summing along the
        #     Δp_l dimension.
        #     Note that the paper implies a sum along Δp_l terms as in
        #     eq. 29, but doesn't actually show them.
        for ph_name, ph_val in signal_pathway.items():
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
                    correl,
                    ax=ax_list[2],
                    human_units=False,
                )
                ax_list[2].set_title(
                    "correlation function (v)\n(includes mask and sum along"
                    " $\\Delta p_l$)"
                )
        # Find optimal f shift based on max of correlation function
        if max_shift is not None:
            delta_f_shift = (
                correl[direct:(-max_shift, max_shift)]
                .run(np.real)
                .argmax(direct)
            )
        else:
            delta_f_shift = correl.run(np.real).argmax(direct)
        # Take s_jk, which is the raw data that has potentially be
        # smooshed, and apply the shift
        s_jk *= np.exp(-1j * 2 * np.pi * delta_f_shift * s_jk.fromaxis(direct))
        f_shift += (
            delta_f_shift  # accumulate all the shifts applied to s_jk to date
        )
        # TODO ☐: this is incorrect, because the mask has not been
        #         applied! (this is not a problem w/ AG changes -- it's
        #         pre-existing)
        #         (Note that once we have masking
        #         functions, we would mod square our
        #         data, and then apply the mask to
        #         the mod squared data -- because
        #         the mask is always real, this is
        #         equivalent to multiplying one
        #         copy of the function by the
        #         mask, then multiplying by an
        #         unmasked copy of the data in
        #         order to calculate our masked
        #         square)
        s_aligned = s_jk.C  # it's probably cheaper to make a copy than to ift
        s_aligned.ft(direct)
        if fl and my_iter == 0:
            fl.image(
                s_aligned,
                ax=ax_list[3],
                human_units=False,
            )
            ax_list[3].set_title("after correlation")
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
            abs(s_aligned.C.sum("repeats")) ** 2
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
        fl.pop_marker()
    # {{{ Make sure returned f_shift is the same ndshape that the
    #     data was originally fed in - meaning we must chunk it back
    #     into the original repeat_dims or rename back to the original
    #     indirect name
    if len(repeats) > 1:
        f_shift.chunk("repeats", repeats)
    else:
        f_shift.rename("repeats", repeats[0])
    # }}}
    return f_shift, sigma, this_mask
