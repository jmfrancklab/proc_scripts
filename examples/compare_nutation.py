import pyspecdata as psd
import pyspecProcScripts as pscr
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def clock_correct(s, axis_along, direct="t2", max_cyc=0.5):
    for_correct = s.C
    Delta = np.diff(s[axis_along][np.r_[0, -1]]).item()
    correction_axis = psd.nddata(np.r_[-0.5:0.5:300j] * max_cyc / Delta, "correction")
    for_correct = for_correct * np.exp(
        -1j * 2 * np.pi * correction_axis * for_correct.fromaxis(axis_along)
    )
    # {{{ determine the best sign flip for each correction
    for j in range(for_correct.shape["correction"]):
        thesign = pscr.determine_sign(for_correct["correction", j])
        for_correct["correction", j] *= thesign
    # }}}
    for_correct.sum(direct)
    return for_correct.sum(axis_along).run(abs).argmax("correction").item()


signal_range = (-250, 250)
fig1, ax1 = plt.subplots(2, 2)
fig2, ax2 = plt.subplots(2, 2)
fig3, ax3 = plt.subplots(2, 2)
fig4, ax4 = plt.subplots(2, 2)
with psd.figlist_var() as fl:
    for searchstr, exptype, nodename, label, axs in [
        [
            "240805_amp0p05_27mM_TEMPOL_nutation",
            "ODNP_NMR_comp/nutation",
            "nutation_1",
            "amplitude 0.05",
            (0,0)
        ],
        [
            "240805_amp0p1_27mM_TEMPOL_nutation",
            "ODNP_NMR_comp/nutation",
            "nutation_1",
            "amplitude = 0.1",
            (0,1)
        ],
        [
            "240805_amp0p2_27mM_TEMPOL_nutation",
            "ODNP_NMR_comp/nutation",
            "nutation_1",
            "amplitude = 0.2",
            (1,0)
        ],
        [
            "240805_amp1_27mM_TEMPOL_nutation",
            "ODNP_NMR_comp/nutation",
            "nutation_1",
            "amplitude = 1",
            (1,1)
        ],
    ]:
        s = psd.find_file(searchstr, exp_type=exptype, expno=nodename, lookup=pscr.lookup_table)
        signal_pathway = s.get_prop("coherence_pathway")
        fl.next("Raw Freq", fig=fig1)
        fig1.suptitle = "compare Raw DCCT"
        ax1[axs[0]][axs[1]].set_title(label)
        if "nScans" in s.dimlabels:
            s.mean("nScans")
        mysgn = pscr.determine_sign(
            pscr.select_pathway(s["t2":signal_range], s.get_prop("coherence_pathway"))
        )
        # {{{ apply overall zeroth order
        s /= pscr.zeroth_order_ph(
                pscr.select_pathway(s['t2':0],s.get_prop("coherence_pathway"))
                )
        d_raw = s.C
        fl.image(pscr.select_pathway(s['t2':(-500,500)].C, signal_pathway), ax = ax1[axs[0]][axs[1]]) 
        d_unc = pscr.select_pathway(s['t2':signal_range].C.real,s.get_prop('coherence_pathway')).integrate('t2')
        # }}}
        # {{{ force zeroth_order on individual betas before 1st order correction
        for j in range(len(s.getaxis("beta"))):
            ph0 = pscr.zeroth_order_ph(
                pscr.select_pathway(
                    s["beta", j]["t2":signal_range], s.get_prop("coherence_pathway")
                )
            )
            s["beta", j] /= ph0
        fl.next("Sign Flipped", fig=fig2)
        fig2.suptitle = "compare Sign Flipped"
        d0 = pscr.select_pathway(s['t2':signal_range].real.C,s.get_prop('coherence_pathway')).integrate('t2')

        ax2[axs[0]][axs[1]].set_title(label)
        # }}}
        # {{{ define mysgn based on phase difference between d_unc and d0
        mysign = (d0/d_unc).angle/np.pi
        mysign = np.exp(1j*np.pi*mysign.run(np.round))
        d_raw *= mysign
        fl.image(pscr.select_pathway(d_raw, s.get_prop("coherence_pathway")), ax=ax2[axs[0]][axs[1]])
        s = d_raw
        s.ift("t2")
        s["t2"] -= s.getaxis("t2")[0]
        s.setaxis(
            "t2", lambda x: x - s.get_prop("acq_params")["tau_us"] * 1e-6
        ).register_axis({"t2": 0})
        s /= pscr.zeroth_order_ph(pscr.select_pathway(s["t2":0.0], signal_pathway))
        s = s["t2":(0, None)]
        s *= 2
        s["t2", 0] *= 0.5
        s.ft("t2")
        s *= mysign
        fl.next("Phased", fig=fig3)
        fig3.suptitle = "compare Phased"
        fl.image(pscr.select_pathway(s, s.get_prop("coherence_pathway")), ax=ax3[axs[0]][axs[1]])
        ax3[axs[0]][axs[1]].set_title(label)
        if "nScans" in s.dimlabels:
            s.mean("nScans")
        fl.next("Integrated", fig=fig4)
        fig4.suptitle = "compare Integrated"
        fl.next('amplitude 0.05 FID vs SE nutation') 
        fin = pscr.select_pathway(s["t2":signal_range].real, signal_pathway)
        print(psd.ndshape(fin))
        fin.integrate("t2"),
        fin.set_error(None)
        A, R, beta_ninety, beta = sp.symbols("A R beta_ninety beta", real = True)
        fl.plot(
            fin,'o',
            color =color,
            ax = ax4[axs[0]][axs[1]]
        )
        ax4[axs[0]][axs[1]].set_title(label)
        f = psd.lmfitdata(fin)
        this_b = 2.5e-5
        this_r = 4e3
        if ("0.1" in label) or ("0.2" in label):
            this_b = 2e-5
            this_r = 3e3
        elif "FID" in nodename:
            this_b = 2.5e-5
            this_r = 3e3
        f.functional_form = functional_form
        f.set_guess(
                A = dict(value=fin.data.max()*1.2, min = fin.data.max()*0.1, max = fin.data.max()*1.2),
                beta_ninety = dict(value=this_b, min = 0, max = 1e3),
                R = dict(value = this_r, min = 0, max = 3e4) 
                )
        f.fit()
        fit = f.eval()
        fl.plot(fit, color = color, ax = ax4[axs[0]][axs[1]])
        psd.gridandtick(plt.gca())
        plt.axvline(22.4e-6)
        ax4[axs[0]][axs[1]].grid()
