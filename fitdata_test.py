from pyspecdata import *
from proc_scripts import hermitian_function_test,zeroth_order_ph
init_logging(level='debug')

#{{{ code bit from proc_IR
from scipy.optimize import leastsq,minimize
#from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
#}}}


#{{{ code for import, pre-fit proc
# {{{ input parameters
date = '200212'
id_string = 'IR_3_30dBm'
clock_correction = 1.785
nodename = 'signal'
filter_bandwidth = 5e3
coh_sel = {'ph1':0,
        'ph2':1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
# }}}
filename = date+'_'+id_string+'.h5'
s = nddata_hdf5(filename+'/'+nodename,
        directory = getDATADIR(exp_type = 'test_equip' ))
s.reorder(['ph2','ph1']).set_units('t2','s')
s *= exp(-1j*s.fromaxis('vd')*clock_correction)
s.ft(['ph2','ph1'])
s.ft('t2', shift=True)
s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
s.ift('t2')
rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
s.setaxis(t2-rough_center)
#s = s['t2':(-25e-3,25e-3)] # select only 50 ms in time domain, because it's so noisy
residual,best_shift = hermitian_function_test(s[
    'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
print("best shift is",best_shift)
s.ft('t2')
s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
s.ift('t2')
ph0 = s['t2':0]['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
print(ndshape(ph0))
if len(ph0.dimlabels) > 0:
    assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
    ph0 = zeroth_order_ph(ph0, fl=None)
    print('phasing dimension as one')
else:
    print("there is only one dimension left -- standard 1D zeroth order phasing")
    ph0 = ph0/abs(ph0)
s /= ph0
s.ft('t2')
s.convolve('t2',10)
s.ift('t2')
s = s['t2':(0,None)]
s *= -1
s['t2',0] *= 0.5
s.ft('t2')
fl.next('signal vs. vd')
s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
s_sliced.sum('t2')
fl.plot(s_sliced,'o')
#}}}
print(ndshape(s_sliced))
print("BEGINNING T1 CURVE...")
f = fitdata(s_sliced)
M0,Mi,T1,vd = sympy.symbols("M_0 M_inf T_1 vd", real=True)
f.functional_form = Mi + (M0-Mi)*sympy.exp(-vd/T1)
print("Functional form",f.functional_form)
print("Function string",f.function_string)
f.fit_coeff = r_[-1,1,1]
#print("f is ",lsafen(f))
fl.next('t1 test')
fl.plot(f,'o',label=f.name())
print("symbolic variables:",f.symbolic_vars) # so I know the ordering for the next step
f.fit()
fl.plot(f.eval(100),label='%s fit'%f.name())
text(0.75, 0.25, f.latex(), transform=gca().transAxes, size='large',
        horizontalalignment='center', color = 'k')
print("output:",f.output())
print("latex:",f.latex())
fl.show();quit()
