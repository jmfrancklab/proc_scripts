import pyspecdata as psd
import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
    frq_mask_fn,
    Delta_p_mask_fn,
    tol=1e-4,
    repeat_dims=[],
    non_repeat_dims=[],
    fig_title="correlation alignment",
    signal_pathway=None,
    max_shift=100.0,
    direct="t2",
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
        Data given in the frequency (not time) domain and the
        phase-cycling (not coherence transfer) domain.
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
    frq_mask_fn : func
        A function which takes nddata and returns a copy that has been
        multiplied by the square root of the frequency-domain mask (see
        DCCT paper).
    Delta_p_mask_fn : func
        A function which takes the 3D data which we call leftbracket
        (:math:`s_{m,n}` in the DCCT paper), and applies the mask over
        the :math:`\\Delta p` (coherence transfer) dimension,
        as well as a sum over :math:`\\Delta p`.
    fl : boolean
        fl=fl to show the plots and figures produced by this
        function otherwise, fl=None.

    Returns
    =======
    f_shift:    np.array
                The optimized frequency shifts for each transient which will
                maximize their correlation amongst each other, thereby aligning
                them.
    """
    # {{{ explicitly check for old arguments, and be generous with format of
    #     input arguments.  Then throw errors for stuff we don't like.
    assert indirect_dim is None, (
        "We updated the correlation function to no longer take indirect_dim as"
        " a kwarg! Now indirect_dim roughly corresponds to repeat_dims"
    )
    assert avg_dim is None, (
        "We updated the correlation function to no longer take avg_dim as a"
        " kwarg!"
    )
    if signal_pathway is None:
        signal_pathway = s_orig.get_prop("coherence_pathway")
    assert signal_pathway is not None, (
        "You need to tell me what the signal pathway is since your data"
        " doesn't have this property set - this is a problem!!"
    )
    if isinstance(non_repeat_dims, str):
        non_repeat_dims = [non_repeat_dims]
    if isinstance(repeat_dims, str):
        repeat_dims = [repeat_dims]
    assert (
        type(repeat_dims) is list and len(repeat_dims) > 0
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
    # }}}
    phcycdims = [j for j in s_orig.dimlabels if j.startswith("ph")]
    if ("nScans" in s_orig.dimlabels) and ("nScans" not in repeat_dims):
        repeat_dims.append("nScans")
    # Make sure the ordering of the repeat_dims matches their order in the
    # original data (this makes smoosh and chunk easier)
    repeat_dims = [j for j in s_orig.dimlabels if j in repeat_dims]
    # {{{ Check there are no left over dimensions unaccounted for by the
    #     direct, phase cycling and declared repeat dimensions
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
    # }}}
    # s_jk below ends up with three dimensions (j = align_dim, k = phcyc and
    # direct nu) and is NOT conj
    if len(repeat_dims) > 1:
        # If there is more than one repeat dim, smoosh into one dimension
        s_jk = s_orig.C.smoosh(repeat_dims, "repeats")
    else:
        s_jk = s_orig.C.rename(
            repeat_dims[0], "repeats"
        )  # even if there isn't an indirect to smoosh we will later be
        #    applying modifications to s_jk that we don't want applied to
        #    s_orig
    s_jk.reorder([direct], first=False)
    for phnames in signal_pathway.keys():
        assert not s_orig.get_ft_prop(phnames), (
            str(phnames) + " must not be in the coherence domain"
        )
    assert s_orig.get_ft_prop(
        direct
    ), "direct dimension must be in the frequency domain"
    ph_len = {j: psd.ndshape(s_orig)[j] for j in signal_pathway.keys()}
    N = s_jk.shape["repeats"]
    # {{{ s_jk is preserved without the mask, and accepts the
    #     progressive shifts to become aligned
    #     s_leftbracket is called that because it becomes (below)
    #     the left square brackets of eq. 29. in Beaton 2022
    #     At this stage, s_mn is equal to s_jk.
    s_leftbracket = frq_mask_fn(s_jk)
    sig_energy = (abs(s_leftbracket) ** 2).data.sum().item() / N
    # }}}
    if fl:
        fl.push_marker()
        fig = fl.next("Correlation Diagnostics")
        fig.suptitle(
            " ".join(
                ["Correlation Diagnostic"]
                + [j for j in [fl.basename] if j is not None]
            )
        )
        gs = GridSpec(1, 4, figure=fig, left=0.05, right=0.95)
        psd.DCCT(
            s_jk,
            fig,
            title="Before correlation\n sig. energy=%g" % sig_energy,
            bbox=gs[0, 0],
        )
    energy_diff = 1.0
    energy_vals = []
    # E_of_avg is the energy calculated from the averaged signal
    # (vs. sig_energy above, which is the energy of the *un*averaged
    # signal.
    # We divide by an extra N because if the signals along the repeats
    # are the same, then the energy of the resulting sum should increase by N
    # (vs taking the square and summing which is what we do for calculating the
    # sig_energy above)
    # TODO ☐ (later): this does not select the coherence pathway, which we want
    #         to do. HOWEVER -- don't do that right now.  Make sure that
    #         everything else works, leaving this todo in until the very
    #         end.  Then, make a comment on github noting that this is
    #         the only remaining todo.  (This can be solved by accepting
    #         data that's in the frequency domain and the coherence
    #         transfer domain, and then just ft/ift'ing both dimensions
    #         together)
    E_of_avg = (
        abs(s_leftbracket.C.sum("repeats")) ** 2
    ).data.sum().item() / N**2
    energy_vals.append(E_of_avg / sig_energy)
    last_E = None
    assert s_jk.get_ft_prop(
        direct
    ), "direct dimension must be in the frequency domain"
    assert not s_jk.get_ft_prop(
        list(s_jk.get_prop("coherence_pathway"))[0]
    ), "phase cycling dimension must not be in the coherence transfer domain!"
    f_shift = 0
    for my_iter in range(100):
        # Note that both s_jk and s_leftbracket
        # change every iteration, because the
        # *data* is updated with every iteration
        logging.debug(psd.strm("*** *** ***"))
        logging.debug(
            psd.strm("CORRELATION ALIGNMENT ITERATION NO. ", my_iter)
        )
        logging.debug(psd.strm("*** *** ***"))
        # note that the frequency mask is applied either (for the first
        # iteration) in the code above or (for subsequent iterations) at
        # the bottom of the for loop
        # {{{ move both the unmasked and masked data into the time domain
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
        # {{{ this applies the Fourier transform from Δφ to Δpₗ
        #     that is found inside the left square bracket of eq. 29.
        #     The paper implies a sum along Δpₗ terms as in eq. 28, but
        #     doesn't actually show them. (Simultaneously maximize all
        #     coherence pathways that are not masked)
        #     Note that because only the left square bracket depends on
        #     Δpₗ, we can apply the coherence mask here, before
        #     multiplication, in order to decrease the dimensionality of
        #     the correlation function.
        # The conjugation needs to be taken before the FT (while in the phase
        # cycling domain so we need to conjugate here (eq 29 in DCCT paper -
        # Beaton 2022)
        s_leftbracket.run(np.conj)
        for ph_name, ph_val in signal_pathway.items():
            s_leftbracket.ft(["Delta%s" % ph_name.capitalize()])
        s_leftbracket = Delta_p_mask_fn(s_leftbracket)
        # }}}
        # the sum over m in eq. 29 only applies to the left bracket,
        # so we just do it here
        correl = s_leftbracket.mean("repeats") * s_jk
        correl.reorder(["repeats", direct], first=False)
        if my_iter == 0:
            logging.debug(psd.strm("holder"))
            if fl:
                correl.reorder([direct], first=False)
                psd.DCCT(
                    correl,
                    fig,
                    title="Correlation function(t)\n (includes ν mask)",
                    bbox=gs[0, 1],
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
                psd.DCCT(
                    correl,
                    fig,
                    title=(
                        "Correlation function (v)\n(includes mask and"
                        " sum along $\\Delta p_l$)"
                    ),
                    bbox=gs[0, 2],
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
        # Take s_jk, which is data that unmasked, but which has all of the
        # shifts from all the previous iterations applied, and apply the
        # shift for this iteration
        s_jk *= np.exp(-1j * 2 * np.pi * delta_f_shift * s_jk.fromaxis(direct))
        # we need to accumulate the total shift
        f_shift += delta_f_shift
        # move back into the frequency domain to analyze the result
        s_jk.ft(direct)
        # the frequency-masked signal (called s_leftbracket here) is not
        # only used to calculate the energy at the end of the for block
        # here, but is also used once we return to the start of the
        # block
        s_leftbracket = frq_mask_fn(s_jk)
        if fl and my_iter == 0:
            psd.DCCT(s_jk, fig, title="After First Iteration", bbox=gs[0, 3])
        logging.debug(
            psd.strm(
                "signal energy per transient (recalc to check that it stays"
                " the same):",
                (abs(s_leftbracket**2).data.sum().item() / N),
            )
        )
        # {{{ Calculate energy difference from last shift to see if
        #     there is any further gain to keep reiterating
        E_of_avg = (
            abs(s_leftbracket.C.sum("repeats")) ** 2
        ).data.sum().item() / N**2
        energy_vals.append(E_of_avg / sig_energy)
        logging.debug(
            psd.strm("averaged signal energy (per transient):", E_of_avg)
        )
        if last_E is not None:
            energy_diff = (E_of_avg - last_E) / sig_energy
            logging.debug(psd.strm(energy_diff))
            if abs(energy_diff) < tol and my_iter > 4:
                break
        # }}}
        last_E = E_of_avg
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
    return f_shift
