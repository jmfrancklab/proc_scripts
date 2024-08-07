from pylab import *
from pyspecdata import *
from scipy.optimize import minimize
from pyspecProcScripts import *
from pyspecProcScripts import fid_from_echo, lookup_table
from sympy import symbols, Symbol, latex
import sympy as sp
from numpy import *
import matplotlib as mpl
def clock_correct(s, axis_along, direct="t2", max_cyc=0.5):
    for_correct = s.C
    Delta = np.diff(s[axis_along][np.r_[0, -1]]).item()
    correction_axis = nddata(np.r_[-0.5:0.5:300j] * max_cyc / Delta, "correction")
    for_correct = for_correct * np.exp(
        -1j * 2 * np.pi * correction_axis * for_correct.fromaxis(axis_along)
    )
    # {{{ determine the best sign flip for each correction
    for j in range(for_correct.shape["correction"]):
        thesign = determine_sign(for_correct["correction", j])
        for_correct["correction", j] *= thesign
    # }}}
    for_correct.sum(direct)
    return for_correct.sum(axis_along).run(abs).argmax("correction").item()
signal_range = (-250,250)
with figlist_var() as fl:
    for searchstr,exptype,nodename, fl.basename  in [
        ['240805_amp0p05_27mM_TEMPOL_nutation','ODNP_NMR_comp/nutation','nutation_1', 'SE amplitude 0.05'],
        ['240805_amp0p1_27mM_TEMPOL_nutation','ODNP_NMR_comp/nutation','nutation_1', 'SE amplitude = 0.1'],
        ['240805_amp0p2_27mM_TEMPOL_nutation','ODNP_NMR_comp/nutation','nutation_1', 'SE amplitude = 0.2'],
        ['240805_amp1_27mM_TEMPOL_nutation','ODNP_NMR_comp/nutation','nutation_1', 'SE amplitude = 1'],
        ]:
        s = find_file(searchstr,exp_type=exptype,expno=nodename, lookup = lookup_table)
        signal_pathway = s.get_prop("coherence_pathway")
        fig, ax_list = plt.subplots(3,1)
        fl.next('Raw Freq',fig = fig)
        fig.suptitle = fl.basename
        if 'nScans' in s.dimlabels:
            fl.image(select_pathway(s.C.mean('nScans'),signal_pathway),ax = ax_list[0])
        else:
            fl.image(select_pathway(s.C,signal_pathway),ax=ax_list[0])
        ax_list[0].set_title('Raw')    
        s.ift('t2')
        s.ft('t2')
        if 'nScans' in s.dimlabels:
            s.mean('nScans')
        mysgn = determine_sign(
                select_pathway(s["t2":signal_range],
                    s.get_prop('coherence_pathway'))
                )
        # {{{ clock correction
        total_corr = 0
        for j in range(5):
            corr = clock_correct(select_pathway(s,s.get_prop('coherence_pathway'))
                    * mysgn
                    *np.exp(-1j*2*pi*total_corr*s.fromaxis('beta')),
                    "beta"
                    )
            total_corr += corr
        s *= np.exp(-1j*2*pi*total_corr*s.fromaxis('beta'))
        mysigns = determine_sign(
                select_pathway(s["t2":signal_range],
                    s.get_prop('coherence_pathway'))
                )
        # {{{ force zeroth_order on individual betas before 1st order correction
        for j in range(len(s.getaxis('beta'))):
            ph0 = zeroth_order_ph(
                    select_pathway(
                        s['beta',j]['t2':signal_range], s.get_prop('coherence_pathway')
                        )
                    )
            s['beta',j] /= ph0
        fl.image(select_pathway(s,s.get_prop('coherence_pathway')), ax = ax_list[1])
        ax_list[1].set_title('same sign')
        #}}}
        s.ift('t2')
        s['t2'] -= s.getaxis('t2')[0]
        s.setaxis('t2', lambda x: x - s.get_prop('acq_params')['tau_us']*1e-6).register_axis({'t2':0})
        s /= zeroth_order_ph(select_pathway(s['t2':0.0],signal_pathway))
        s = s['t2':(0,None)]
        s *= 2
        s['t2',0] *= 0.5
        s.ft('t2')
        s *= -mysigns
        fl.image(select_pathway(s, s.get_prop('coherence_pathway')),ax = ax_list[2])
        ax_list[2].set_title('Phased and flipped')
        fl.basename = None
        fl.next('integrated')
        if 'nScans' in s.dimlabels:
            s.mean('nScans')
        fl.plot(select_pathway(s['t2':signal_range],signal_pathway).C.real.integrate('t2'),
                'o',label = fl.basename, human_units = False)


