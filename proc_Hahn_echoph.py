from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from sympy import symbols
fl = fl_mod()
t2 = symbols('t2')
logger = init_logging('info')
def select_pathway(s,pathway):
    retval = s
    for k, v in pathway.items():
        retval = retval[k,v]
    return retval
signal_pathway = {'ph1':1,'ph2':0}
for searchstr, exp_type, nodename, postproc, label_str, slice_f in [
        ('200302_alex_probe_water', 'test_equip', 'signal', 
            'spincore_Hahn_echoph_v1','microwaves off',(-2.5e3,2.5e3)),
        ]:
    #{{{loads raw data and plots
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict,fl=fl)
    s.mean('nScans')    
    #}}}
    #{{{rough centering of sliced data 
    s = s['t2':slice_f]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    logger.debug(strm(ndshape(s)))
    #}}}
    #{{{ apply phase corrections
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))
    logger.info(strm("best shift is",best_shift))
    s_uncorrected = s.C.ft('t2')
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0},nearest=False)
    ph0 = s['t2':0]['ph2',0]['ph1',1]
    logger.info(strm(ndshape(ph0)))
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0, fl=fl)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = ph0/abs(ph0)
    s /= ph0
    #}}}
    #{{{visualizes the data after hermitian function test and phasing 
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2', pad=512) # power of 2 FT
    fl.image(s.C.convolve('t2',10)) # so that resolution of
    # plot isn't higher than that of screen -- note that we do
    # this on a copy of the data, since we don't actually want
    # to alter the data here
    #}}}
    #{{{slice out FID from echo
    s.ift('t2')
    s = slice_FID_from_echo(s)
    s.ft('t2')
    #}}}
    #{{{visualize final processed data
    s = select_pathway(s,signal_pathway)
    fl.next('processed data')
    fl.plot(s_uncorrected['ph2',-2]['ph1',1],'ko',
            label='without time-axis correction')
    fl.plot(s,'ro',label='with time-axis correction')
    fl.next('Spectrum FT')
    fl.plot(s.real, alpha=0.5, label='real - %s'%label_str)
    fl.plot(s.imag, alpha=0.5, label='imag - %s'%label_str)
    #}}}
fl.show()
