from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
logger = init_logging("info")
# {{{ input parameters
clock_correction = 1.785
filter_bandwidth = 5e3
coh_sel = {'ph1':0,
        'ph2':1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
# }}}

for searchstr,exp_type,nodename, postproc in [
        ('200212_IR_3_30dBm', 'test_equip', 'signal', 
            'spincore_IR_v1'),
        ]:
    fl.basename = searchstr
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename,
            postproc=postproc, lookup=postproc_dict,
            dimname='indirect')
    s *= exp(-1j*s.fromaxis('indirect')*clock_correction)
    #logger.info(strm(s.dimlabels))
    #{{{rough centers data
    #fl.next('filtered + rough centered data')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    #}}}
    #{{{hermitian function test and apply best shift
    fl.next('frequency domain before')
    fl.image(s)
    s.ift('t2')
    best_shift = hermitian_function_test(s[
        'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
    logger.info(strm("best shift is",best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    s.reorder(['ph2','ph1','indirect'])
    fl.next('time domain after hermitian test')
    fl.image(s)
    fl.next('frequency domain after')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    #}}}
    #{{{zeroth order phase correction
    ph0 = s['t2':0]['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
    logger.info(strm(ndshape(ph0)))
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0, fl=fl)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.image(s.C.convolve('t2',10))
    #}}}
    #{{{select t2 axis range and 
    s.ift('t2')
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    #}}}
    #{{{decay curve and fitting
    fl.next('signal vs. vd')
    s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]*-1 # bc inverted at high powers
    s_sliced.sum('t2')
    fl.plot(s_sliced,'o')
    logger.info(strm(ndshape(s_sliced)))
    logger.info(strm("BEGINNING T1 CURVE..."))
    s = fitdata(s_sliced)
    M0,Mi,R1,vd = sympy.symbols("M_0 M_inf R_1 indirect", real=True)
    s.functional_form = Mi + (M0-Mi)*sympy.exp(-vd*R1)
    logger.info(strm("Functional form", s.functional_form))
    # JF notes that we want to be able to set the guess using a dictionary
    # here (which is what I think setting fit_coeff was doing), then plot
    # the guess to make sure that's what we're doing -- like so
    fl.next('t1 test')
    s.set_guess(Mi=-1, M0=1, R1=1)
    # work, currently -- we will need a pull request on pyspecdata as well
    # to make it work
    fl.plot(g, 'o', label="data")
    s.settoguess()
    fl.plot(s, '-', label='initial guess')
    s.fit()
    fl.plot(s.eval(100),label='fit')
    text(0.75, 0.25, s.latex(), transform=gca().transAxes, size='large',
            horizontalalignment='center',color='k')
    print("output:",s.output())
    print("latex:",s.latex())
fl.show()

