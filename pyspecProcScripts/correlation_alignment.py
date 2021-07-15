from pylab import*
from pyspecdata import *
from scipy.optimize import leastsq, minimize, basinhopping
import numpy as np
from matplotlib.ticker import FuncFormatter
import logging
@FuncFormatter
def to_percent(y, position):
    s = '%.2g'%(100 * y)
    if rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'
def correl_align(s, align_phases=False,tol=1e-4,indirect_dim='indirect',
        fig_title='correlation alignment',signal_pathway = {'ph1':0,'ph2':1}, 
        shift_bounds=False, avg_dim = None, max_shift = 100., sigma=20.,direct='t2',fl=None):
    """
    Align transients collected with chunked phase cycling dimensions along an indirect
    dimension based on maximizing the correlation across all the transients and repeat
    alignment until the calculated signal energy remains constant to within a given 
    tolerance level.

    Parameters
    ==========
    s:  nddata
        A nddata object which contains phase cycle dimensions and an
        indirect dimension.
    align_phases:   boolean
    tol:            float
                    Sets the tolerance limit for the alignment procedure.
    indirect_dim:   str
                    Name of the indirect dimension along which you seek to align
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
                    Sigma value for the Gaussian fitting. Related to the line width
                    of the given data.
    fl:             boolean 
                    fl=fl to show the plots and figures produced by this function
                    otherwise, fl=None.
                    
    Returns
    =======
    f_shift:    array
                The optimized frequency shifts for each transient which will 
                maximize their correlation amongst each other, thereby aligning
                them.
    sigma:      float
                The width of the Gaussian function used to frequency filter
                the data in the calculation of the correlation function.
    """

    logging.info(strm("Applying the correlation routine"))
    if avg_dim:
        phcycdims = [j for j in s.dimlabels if j.startswith('ph')]
        indirect = set(s.dimlabels)-set(phcycdims)-set([direct])
        indirect = [j for j in s.dimlabels if j in indirect]
        avg_dim_len = len(s.getaxis(avg_dim))
        s.smoosh(indirect)
    s.ift(list(signal_pathway.keys()))
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    ph_len = {j:ndshape(s)[j] for j in signal_pathway.keys()}
    N = ndshape(s)[indirect_dim]
    sig_energy = (abs(s)**2).data.sum().item() / N
    if fl:
        fig_forlist, ax_list = plt.subplots(1, 5, figsize=(7,7))
        fl.next("Correlation Diagnostics")
        fig_forlist.suptitle(" ".join(["Correlation Diagnostic"] + [j for j in [fl.basename] if j is not None]))
        fl.image(s.C.setaxis(indirect_dim,'#').set_units(indirect_dim,'scan #'),ax=ax_list[0],human_units=False)
        ax_list[0].set_title('before correlation\nsig. energy=%g'%sig_energy)
    energy_diff = 1.
    i = 0
    energy_vals = []
    this_E = (abs(s.C.sum(indirect_dim))**2).data.sum().item() / N**2
    energy_vals.append(this_E / sig_energy)
    last_E = None
    for_nu_center =s.C
    for_nu_center.ft(list(signal_pathway))
    for x in range(len(signal_keys)):
        for_nu_center = for_nu_center[signal_keys[x],signal_values[x]]
    nu_center = for_nu_center.mean(indirect_dim).C.argmax('t2')
    logging.info(strm("Center frequency", nu_center))
    for my_iter in range(100):
        i += 1
        logging.info(strm("*** *** ***"))
        logging.info(strm("CORRELATION ALIGNMENT ITERATION NO. ",i))
        logging.info(strm("*** *** ***"))
        if align_phases:
            ph0 = s.C.sum('t2')
            ph0 /= abs(ph0)
            s /= ph0
        s.ift('t2')
        s_copy = s.C
        s_copy.ft('t2')
        s_copy *= exp(-(s_copy.fromaxis('t2')-nu_center)**2/(2*sigma**2))
        s_copy.ift('t2')
        s_copy2 = s.C
        for k,v in ph_len.items():
            ph = ones(v)
            s_copy *= nddata(ph,'Delta'+k.capitalize())
            s_copy.setaxis('Delta'+k.capitalize(),'#')
        correl = s_copy * 0 
        for k,v in ph_len.items():
            for ph_index in range(v):
                s_copy['Delta%s'%k.capitalize(),ph_index] = s_copy['Delta%s'%k.capitalize(),
                        ph_index].run(lambda x, axis=None: roll(x,ph_index,axis=axis),k)
        for j in range(1,N):
            correl += s_copy2 * s_copy.C.run(lambda x, axis=None: roll(x,j,axis=axis),
                indirect_dim).run(conj)
        correl.reorder([indirect_dim,'t2'],first=False)
        if my_iter ==0:
            logging.info(strm("holder"))
            if fl:
                fl.image(correl.C.setaxis(indirect_dim,'#').set_units(indirect_dim,'scan #'),
                        ax=ax_list[1])
                ax_list[1].set_title('correlation function (t), \nafter apod')
        correl.ft_clear_startpoints('t2')
        correl.ft('t2', shift=True, pad=2**14)
        for k,v in signal_pathway.items():
            correl.ft(['Delta%s'%k.capitalize()])
            correl = correl['Delta'+k.capitalize(),v]+correl['Delta'+k.capitalize(),0]
        if my_iter ==0:
            logging.info(strm("holder"))
            if fl:
                fl.image(correl.C.setaxis(indirect_dim,'#').set_units(indirect_dim,'scan #'),
                        ax=ax_list[2],human_units=False)
                ax_list[2].set_title('correlation function (v), \nafter apod')
        if shift_bounds:
            f_shift = correl['t2':(-max_shift,max_shift)].run(real).argmax('t2')
        else:
            f_shift = correl.run(real).argmax('t2')
        s_copy = s.C
        s_copy *= exp(-1j*2*pi*f_shift*s_copy.fromaxis('t2'))
        s.ft('t2')
        s_copy.ft('t2')
        if my_iter ==0:
            logging.info(strm("holder"))
            if fl:
                fl.image(s_copy.C.setaxis(indirect_dim,'#').set_units(indirect_dim,'scan #'),
                        ax=ax_list[3],human_units=False)
                ax_list[3].set_title('after correlation\nbefore ph0 restore')
        logging.info(strm('signal energy per transient (recalc to check that it stays the same):',(abs(s_copy**2).data.sum().item() / N)))

        this_E = (abs(s_copy.C.sum(indirect_dim))**2).data.sum().item() / N**2
        energy_vals.append(this_E / sig_energy)
        logging.info(strm('averaged signal energy (per transient):', this_E))
        if last_E is not None:
            energy_diff = (this_E - last_E)/sig_energy
            logging.info(strm(energy_diff))
            if abs(energy_diff) < tol and my_iter > 4:
                break
        last_E = this_E
    if fl is not None: 
        fl.next('correlation convergence')
        fl.plot(array(energy_vals),'x')
        gca().yaxis.set_major_formatter(to_percent)
    if fl is not None:
        fl.image(s_copy.C.setaxis(indirect_dim,'#').set_units(indirect_dim,'scan #'),ax=ax_list[4])
        ax_list[4].set_title('after correlation\nph0 restored \nsig. energy=%g'%sig_energy)
    if avg_dim:
        s.chunk(avg_dim,[avg_dim,'power'],[avg_dim_len,-1])
        s.reorder(['ph1',avg_dim,'power','t2'])
    return f_shift,sigma    

