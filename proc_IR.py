from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')

# {{{ input parameters
clock_correction = 1.785
filter_bandwidth = 5e3
coh_sel = {'ph1':0,
        'ph2':1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
# }}}

for searchstr,exp_type,nodename, postproc in [
        ('200303_IR_AER_6','test_equip','signal','spincore_ONDP_v1'),
        ]:
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename,
            postproc=None, lookup=postproc_dict,
            dimname='indirect')
    logger.info(strm(s.dimlabels))
    
    #rough centers data
    fl.next('filtered + rough centered data')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but(
            't2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    
    #hermitian function test and apply best shift
    best_shift = hermitian_function_test(s[
        'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
    logger.info(strm("best shift is",best_shift))
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.reorder(['ph2','ph1','vd'])
    s.ift('t2')
    fl.next('time domain after hermitian test')
    fl.image(s)

    #zeroth order phase correction
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
    s.convolve('t2',10)
    fl.image(s)
    
    #select t2 axis range and 
    s.ift('t2')
    s = s['t2':(0,None)]
    s *= -1
    s['t2',0] *= 0.5
    s.ft('t2')
    
    
    fl.next('signal vs. vd')
    s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]*-1 # bc inverted at high powers
    s_sliced.sum('t2')
    fl.plot(s_sliced,'o')
    print(ndshape(s_sliced))
    print("BEGINNING T1 CURVE...")
    s = fitdata(s_sliced)
    M0,Mi,T1,vd = sympy.symbols("M_0 M_inf T_1 vd", real=True)
    s.functional_form = Mi + (M0-Mi)*sympy.exp(-vd/T1)
    print("Functional form", s.functional_form)
    print("Function string"),s.function_string)
    s.fit_coeff = r_[-1,1,1]
    fl.next('t1 test')
    fl.plot(s, 'o', label=s.name())
    print("symbolic variables:",f.symbolic_vars)
    fl.plot(f.fitfunc_multiarg(-3000,3000,1,s.fromaxis('vd')))
    fl.plot(s.getaxis('vd'),s.fitfunc_raw(r_[-3000,3000,1],s.getaxis('vd')),'--')
    s.fit()
    fl.plot(s.eval(100),label='%s fit'%s.name())
    text(0.75, 0.25, s.latex(), transform=gca().transAxes, size='large',
            horizontalalignment='center',color='k')
    print("output:",s.output())
    print("latex:",s.latex())
    fl.show();quit()

    s_forerror = s['ph2',coh_err['ph2']]['ph1',0]
    # variance along t2 gives error for the mean, then average it across all other dimension, then sqrt for stdev
    s_forerror.run(lambda x: abs(x)**2).mean_all_but(['vd']).run(sqrt)
    s_sliced.mean('t2').set_error(s_forerror.data)
    fl.plot(s_sliced,'o')
    fl.plot(s_sliced.imag,'o')
    fl.next('Spectrum - freq domain')
    fl.plot(s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
    fitfunc = lambda p, x: p[0]*(1-2*exp(-x*p[1]))
    x = s_sliced.fromaxis('vd')
    errfunc = lambda p: fitfunc(p,x).data - s_sliced.data.real
    p_ini = [1.0,10.0]
    p_opt,success = leastsq(errfunc, p_ini[:])
    assert success in [1,2,3], "Fit did not succeed"
    T1 = 1./p_opt[1]
    logger.info(strm("T1:",T1,"s"))
    fl.next('fit')
    fl.plot(s_sliced, 'o', label='data')
    fl.plot(s_sliced.imag, 'o', label='data')
    new_x = nddata(r_[0:s_sliced.getaxis('vd')[-1]:500j],'vd')#.set_units('vd','s')
    fl.plot(fitfunc(p_opt,new_x), label='fit')
    gridandtick(gca())
    fl.show();quit()
