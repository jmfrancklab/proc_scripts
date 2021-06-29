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
        an nddata object which contains phase cycle dimensions and an
        indirect dimension
    align_phases:   boolean
    indirect_dim:   str
                    name of the indirect dimension along which you seek to align
                    the transients
    fig_title:      str
                    name for the figures generated
    tol:            float
                    sets the tolerance limit for the alignment procedure
    signal_pathway: dict
                    dictionary containing the signal pathway
    shift_bounds:   boolean
                    keeps f_shift to be within a specified
                    limit (upper and lower bounds given by max_shift)
                    which should be around the location of the expected
                    signal
    avg_dim:        str
                    dimension along which the data is being averaged
    max_shift:      float
                    specifies the upper and lower bounds to the range over
                    which f_shift will be taken from the correlation function.
                    shift_bounds must be True
    sigma:          int
                    sigma value for the gaussian fitting. Related to the linewidth
                    of the given data.
    fl:             boolean 
                    fl=fl to show the plots and figures produced by this function
                    otherwise, fl=None
                    
    Returns
    =======
    f_shift:    array
                the optimized frequency shifts for each transient which will 
                maximize their correlation amongst each other, thereby aligning
                them
    sigma:      float
                the width of the Gaussian function used to frequency filter
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
    if fl is not None:
        fl.next('before correlation\nsig. energy=%g'%sig_energy + fig_title)
        fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'),human_units=False)
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
        freq_filter = True
        s_copy.ft('t2')
        sigma = sigma
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
        if my_iter == 0:
            logging.info(strm("holder"))
            if fl is not None:
                fl.next('Look at correlation function - time domain')
                fl.image(correl.C.setaxis('vd','#').set_units('vd','scan #'),human_units=False)
        if my_iter ==0:
            logging.info(strm("holder"))
            if fl is not None:
                fl.next('correlation function\ntime domain, after apod')
                fl.image(correl.C.setaxis('vd','#').set_units('vd','scan #'),human_units=False)
        correl.ft_clear_startpoints('t2')
        correl.ft('t2', shift=True, pad=2**14)
        for k,v in ph_len.items():
                correl.ft(['Delta%s'%k.capitalize()])
        for j in range(1,len(ph_len.keys())+1):
            for k,v in ph_len.items():
                correl = correl['Delta'+k.capitalize(),j]+correl['Delta'+k.capitalize(),0]
        if my_iter ==0:
            logging.info(strm("holder"))
            if fl is not None:
                fl.next('correlation function \nfreq domain, after apod')
                fl.image(correl.C.setaxis('vd','#').set_units('vd','scan #'),human_units=False)
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
            if fl is not None:
                fl.next('after correlation\nbefore ph0 restore')
                fl.image(s_copy.C.setaxis('vd','#').set_units('vd','scan #'),human_units=False)
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
        fl.next('after correlation\nph0 restored sig. energy=%g'%sig_energy)
        fl.image(s_copy.C.setaxis('vd','#').set_units('vd','scan #'),human_units=False)
    s.ift('t2')
    s *= np.exp(-1j*2*pi*f_shift*s.fromaxis('t2'))
    s.ft('t2')
    if avg_dim:
        s.chunk(avg_dim,[avg_dim,'power'],[avg_dim_len,-1])
        s.reorder(['ph1',avg_dim,'power','t2'])
    return s,f_shift,sigma    

