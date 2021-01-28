from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import hermitian_function_test, zeroth_order_ph, recovery
from sympy import symbols, latex, Symbol
from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np
from sympy import exp as s_exp
fl = figlist_var()
t2 = symbols('t2')
# {{{ input parameters
date = '210127'
id_string = 'OHTEMPO10mM_cap_probe_IR_1'
nodename = 'signal'
filter_bandwidth = 5e3
coh_sel = {'ph1':0,
        'ph2':1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
# }}}
filename = date+'_'+id_string+'.h5'
s = nddata_hdf5(filename+'/'+nodename,
        directory = getDATADIR(exp_type='test_equip'))
vd_axis = s.getaxis('vd')
s.reorder(['ph2','ph1']).set_units('t2','s')
s.ft('t2',shift=True)
fl.next('raw data')
fl.image(s.C.setaxis('vd','#'))
#s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
fl.next('raw data -- coherence channels')
s.ft(['ph2','ph1'])
fl.image(s.C.setaxis('vd','#'))
s.ift('t2')
fl.next('time domain cropped log')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #').cropped_log())
rough_center = abs(s).convolve('t2',10).mean_all_but('t2').argmax('t2').item()
s.setaxis(t2-rough_center)
fl.next('rough centering')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
s.ft('t2')
s.ift('t2')
best_shift = hermitian_function_test(s[
    'ph2',1]['ph1',0])
logger.info(strm("best shift is", best_shift))
s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
fl.next('time domain after hermitian test')
fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
ph0 = s['t2':0]['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
if len(ph0.dimlabels) > 0:
    assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
    ph0 = zeroth_order_ph(ph0)
    print('phasing dimension as one')
else:
    print("there is only one dimension left -- standard 1D zeroth order phasing")
    ph0 = ph0/abs(ph0)
s /= ph0
fl.next('frequency domain -- after hermitian function test and phasing')
s.ft('t2')
fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
s.ift('t2')
fl.next('check phasing -- real')
fl.plot(s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
gridandtick(plt.gca())
fl.next('check phasing -- imag')
fl.plot(s[
    'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']].imag)
gridandtick(plt.gca())
s = s['t2':(0,None)]
s['t2',0] *= 0.5
fl.next('phased and FID sliced')
fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
fl.next('phased and FID sliced -- frequency domain')
s.ft('t2')
# }}}
fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))
fl.next('recovery curve')
s_signal = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
s_forerror = s['ph2',coh_err['ph2']]['ph1',coh_err['ph1']]
# variance along t2 gives error for the mean, then average it across all other dimension, then sqrt for stdev
print(ndshape(s_forerror))
s_forerror.run(lambda x:abs(x)**2).mean_all_but(['vd','t2']).integrate('t2').run(sqrt)
s_signal.integrate('t2').set_error(s_forerror.data)
fl.plot(s_signal,'o',label='real')
fl.plot(s_sliced.imag,'o',label='imaginary')
fl.next('Spectrum - freq domain')
s = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
fl.plot(s)
#s_sliced *= -1
fitfunc = lambda p, x: p[0]*(1-2*np.exp(-x*p[1]))
x = s_signal.fromaxis('vd')
errfunc = lambda p: fitfunc(p,x).data - s_signal.data.real
p_ini = [1.0,10.0]
p_opt,success = leastsq(errfunc, p_ini[:])
assert success in [1,2,3], "Fit did not succeed"
T1 = 1./p_opt[1]
print("T1:",T1,"s")

f = fitdata(s_signal)
error = fitdata(s_forerror)
M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd", real=True)
error.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
f.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
error.fit()
f.fit()
print("output:",f.output())
print("latex:",f.latex())
T1 = 1./f.output('R_1')
fl.next('fit')
fl.plot(s_signal,'o', label='actual data')
fl.plot(s_signal.imag,'--',label='actual imaginary')
fl.plot(f.eval(100),label='fit')
fl.plot(s_forerror,'o',label='error')
fl.plot(error.eval(100),label='fit error')
plt.text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes,size='large',
        horizontalalignment='center',color='k')
ax = plt.gca()
fl.show();quit()
