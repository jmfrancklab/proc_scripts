from proc_scripts import * 
from pyspecdata import *
from sympy import symbols


def slice_FID_from_echo(s, max_t=None, signal_pathway={'ph1':1, 'ph2':-2}, fl=None):
    """Automatically slices the FID out of an echo.

    Determines the proper position for :math:`t=0` from hermitian function test.
    This is followed by zeroth order phase correcting and slicing from the center of 
    the echo (:math:`t=0` onward to the t_range defined

    Parameters
    ==========
    max_t:  int
            where the slice ends in the time domain.
    ph1:    int
            the selection of ph1 in the coherence pathway
    ph2:    int
            selection of ph2 in coherence pathway

    Returns
    =======
    s:      phase corrected and sliced data
    """
    def select_pathway(s,pathway):
        myslice = s
        for k, v in pathway.items():
            myslice = myslice[k,v]
        return myslice    
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    coh_slice = select_pathway(s['t2':0],signal_pathway)
    print(coh_slice.dimlabels)
    if len(coh_slice.dimlabels) > 0:
        assert len(coh_slice.dimlabels) == 1, repr(ndshape(coh_slice.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(coh_slice, fl=fl)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = coh_slice/abs(coh_slice)
    s /= ph0
    if fl is not None:
        fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range', s, (0,max_t))
    s = s['t2':(0,max_t)]
    s['t2':0] *= 0.5
    if 'power' in s.dimlabels:
        if select_pathway(s['t2':0]['power',0],signal_pathway).real < 0:
            s *= -1
    elif 'vd' in s.dimlabels:
        if s['t2':0]['ph2',ph2]['ph1',ph1]['vd',-1].real < 0:
            s *= -1
    return s
