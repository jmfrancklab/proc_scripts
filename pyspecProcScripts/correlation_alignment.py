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
    repeat_dims=[],
    non_repeat_dims=[],
    fig_title="correlation alignment",
    signal_pathway={"ph1": 0, "ph2": 1},
    max_shift=100.0,
    sigma=20.0,
    direct="t2",
    frq_mask_fn=None,
    ph_mask_fn=None,
    fl=None,
    indirect_dim=None,  # no longer used
    avg_dim=None,  # no longer used
):
    """
    Align transients collected with chunked phase cycling dimensions along an
    indirect dimension based on maximizing the correlation across all the
    transients and repeat alignment until the calculated signal energy remains
    constant to within a given tolerance level.

    This implements Eq 29 in the DCCT paper (Beaton et al JCP 2022),
    remembering that the FT of the cross-correlation is the product
    of the two functions, with the *first* on complex conjugate.

    With both repeat_dims and non_repeat_dims defined, we can use a set
    expression to find the indirect dimension (not phase cycling and not
    direct) and to make sure all dimensions are either
    (1) nScans
    (2) specified as safe to align along (by being listed in `repeat_dims`)
    or (3) listed in the `non_repeat_dims`.

    Parameters
    ==========
    s_orig:  psd.nddata
        A psd.nddata object which contains phase cycle dimensions and an
        indirect dimension.
        This data is not modified.
    tol:            float
                    Sets the tolerance limit for the alignment procedure.
    repeat_dims : list (default [])
        List of the dimensions along which the signal is
        essentially repeated and therefore can be aligned.
        The nScans dimension is automatically added to this list if it exists.
        **Important**: the signal along a dimension of repeats must be the same
        sign.
        Therefore, *e.g.*, if it's derived from inversion recovery data,
        you need to flip the sign →
        That's what :meth:`pyspecProcScripts.phasing.determine_sign` is for.
    nonrepeat_dims : list (default [])
        These are indirect dimensions that we are not OK to align
        (e.g. the indirect dimension of a 2D COSY experiment).
    fig_title : str
        Title for the figures generated.
    signal_pathway : dict
        Dictionary containing the signal pathway.
    max_shift : float
        Specifies the upper and lower bounds to the range over
        which f_shift will be taken from the correlation function.
        Shift_bounds must be True.
        If it's set to None, then no bounds are applied.
    sigma : int
        Sigma value for the Gaussian mask. Related to the line
        width of the given data.
    frq_mask_fn : func
        A function which takes nddata and returns a copy with a frequency
        mask applied that only leaves a bandwidth surrounding the signal as
        nonzero.
    ph_mask_fn : func 
        A function which takes the 3D data which we call leftbracket (and
        pertains to s_{m,n} in the DCCT paper) and filters the selected
        pathwaysover which to sum.
    fl : boolean
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
    # {{{ explicitly check for old arguments
    assert indirect_dim is None, (
        "We updated the correlation function to no longer take indirect_dim as"
        " a kwarg! Now indirect_dim roughly corresponds to repeat_dims"
    )
    assert avg_dim is None, (
        "We updated the correlation function to no longer take avg_dim as a"
        " kwarg!"
    )
    # }}}
    if isinstance(non_repeat_dims, str):
        non_repeat_dims = [non_repeat_dims]
    if isinstance(repeat_dims, str):
        repeat_dims = [repeat_dims]
    assert (
        type(repeat_dims) is list
    ), "the repeat_dims kwarg needs to be a list of strings"
    assert (
        len(repeat_dims) > 0
    ), "You must tell me which dimension contains the repeats!"
    temp = set(repeat_dims) - set(s_orig.dimlabels)
    assert len(temp) == 0, (
        f"{temp} were not found in the data dimensions, but were specified in"
        " `repeat_dims`"
    )
    temp = set(non_repeat_dims) - set(s_orig.dimlabels)
    assert len(temp) == 0, (
        f"{temp} were not found in the data dimensions, but were specified in"
        " `nonrepeat_dims`"
    )
    assert frq_mask_fn is not None, (
        "You need to give me a function that will frequency filter the "
        "signal! It should take the arguments: s, signal_pathway, direct"
        " and indirect"
    )

    phcycdims = [j for j in s_orig.dimlabels if j.startswith("ph")]
    # TODO ✓: check my modifications -- the following was incorrectly called
    # safe_repeat_dims rather than indirect (what it was called before, and a
    # more accurate name) see
    # https://github.com/jmfrancklab/proc_scripts/pull/141#discussion_r2042969536
    if ("nScans" in s_orig.dimlabels) and ("nScans" not in repeat_dims):
        repeat_dims.append("nScans")
    # Make sure the ordering of the repeat_dims matches their order in the
    # original data (this makes smoosh and chunk easier)
    repeat_dims = [j for j in s_orig.dimlabels if j in repeat_dims]
    # Check there are no left over dimensions unaccounted for by the direct,
    # phase cycling and declared repeat dimensions
    temp = (
        set(s_orig.dimlabels)
        - set(phcycdims)
        - set([direct])
        - set(repeat_dims)
        - set(non_repeat_dims)
    )
    assert len(temp) == 0, (
        "Aside from the dimension called nScans, the direct dimension"
        f" {direct}, and dimensions called ph... (for phase cycling), you need"
        " to acccount for all your indirect dimensions!  The following are"
        f" unaccounted for: {temp}"
    )
    # s_jk below ends up with three dimensions (j = align_dim, k = phcyc and
    # direct nu) and is NOT conj
    if len(repeat_dims) > 1:
        # If there is more than one repeat dim, smoosh into one dimension
        s_jk = s_orig.C.smoosh(repeat_dims, "repeats")
    else:
        s_jk = s_orig.C.rename(
            repeat_dims[0], "repeats"
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
    ph_len = {j: psd.ndshape(s_orig)[j] for j in signal_pathway.keys()}
    N = s_jk.shape["repeats"]
    s_jk.ft(list(signal_pathway))
    for_sig_E = frq_mask_fn(s_jk)
    s_jk.ift(list(signal_pathway))
    sig_energy = (abs(for_sig_E) ** 2).data.sum().item() / N
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
        s_jk.ft(list(signal_pathway))
        s_leftbracket = frq_mask_fn(s_jk)
        s_jk.ift(list(signal_pathway))
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
        # {{{ this applies the Fourier transform from Δφ to Δp_l
        #     that is found inside the left square bracket of eq. 29.
        #     Then, it performs a sum that is equivalent to applying a
        #     mask along the Δp_l dimension and then summing along the
        #     Δp_l dimension.
        #     Note that the paper implies a sum along Δp_l terms as in
        #     eq. 29, but doesn't actually show them.
        s_leftbracket = ph_mask_fn(s_leftbracket, signal_pathway)
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
    if len(repeat_dims) > 1:
        f_shift.chunk("repeats", repeat_dims)
    else:
        f_shift.rename("repeats", repeat_dims[0])
    # }}}
    return f_shift, sigma
