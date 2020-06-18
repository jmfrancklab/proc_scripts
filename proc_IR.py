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
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename,
            postproc=postproc, lookup=postproc_dict,
            dimname='indirect')
    s *= exp(-1j*s.fromaxis('indirect')*clock_correction)
    logger.info(strm(s.dimlabels))
    #{{{rough centers data
    fl.next('filtered + rough centered data')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but(
            't2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    #}}}
    #{{{slicing out FID from echo and centering
    s = slice_FID_from_echo(s,1,0) # so, to be clear, this includes all phasing?
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
    M0,Mi,R1,vd = sympy.symbols("M_0 M_inf T_1 indirect", real=True)
    s.functional_form = Mi + (M0-Mi)*sympy.exp(-vd*R1)
    logger.info(strm("Functional form", s.functional_form))
    # JF notes that we want to be able to set the guess using a dictionary
    # here (which is what I think setting fit_coeff was doing), then plot
    # the guess to make sure that's what we're doing -- like so
    fl.next('t1 test')
    s.set_guess({M0:-1, Mi:1, R1:1}) # this is the only line that will not
    # work, currently -- we will need a pull request on pyspecdata as well
    # to make it work
    fl.plot(s, 'o', label=s.name())
    s.settoguess()
    fl.plot(s, '-', label='initial guess')
    s.fit()
    fl.plot(s.eval(100),label='%s fit'%s.name())
    text(0.75, 0.25, s.latex(), transform=gca().transAxes, size='large',
            horizontalalignment='center',color='k')
    print("output:",s.output())
    print("latex:",s.latex())
fl.show()

