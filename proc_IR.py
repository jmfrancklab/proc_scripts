from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from proc_scripts.fitting import recovery
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
    #{{{filter data
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    #}}}
    #{{{hermitian function test and apply best shift
    fl.next('frequency domain before')
    fl.image(s)
    s.ift('t2')
    best_shift = hermitian_function_test(s[
        'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']],fl=fl)
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
    #{{{recovery curve and fitting
    s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]*-1 # bc inverted at high powers
    # below, f_range needs to be defined
    M0,Mi,R1,vd = sympy.symbols("M_0 M_inf R_1 indirect",real=True)
    f,T1,g = recovery(s_sliced, (-100,100),
            guess={M0:-500, Mi:500, R1:1})
    fl.plot_curve(f,'inversion recovery curve',guess=g)
fl.show()

