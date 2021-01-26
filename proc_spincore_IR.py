from pyspecdata import *
from scipy.optimize import leastsq,minimize
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
import matplotlib.pyplot as plt
import numpy as np
fl = figlist_var()
t2 = symbols('t2')
# {{{ input parameters
date = '210120'
id_string = 'OHTEMPO10mM_sol_probe_IR'
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
s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
fl.next('raw data -- coherence channels')
s.ft(['ph2','ph1'])
fl.image(s.C.setaxis('vd','#'))
s.ift('t2')
rough_center = abs(s).convolve('t2',10).mean_all_but('t2').argmax('t2').item()
s.setaxis(t2-rough_center)
s.ft('t2')
fl.next('before clock correction')
fl.image(s.C.setaxis('vd','#'))
clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr') # setaxis not required if passed this way
s_clock = s['ph1',0]['ph2',1].sum('t2')
min_index = abs(s_clock).argmin('vd', raw_index=True).item()
s_clock *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
s_clock['vd',:min_index+1] *= -1
s_clock.sum('vd').run(abs)
fl.next('clock correction')
fl.plot(s_clock,'.',alpha=0.7)
clock_corr = s_clock.argmax('clock_corr').item()
plt.axvline(x=clock_corr, alpha=0.5, color='r')
s *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
fl.next('after auto-clock correction')
fl.image(s.C.setaxis('vd','#'))
#fl.show();quit()
s.ift('t2')
residual,best_shift = hermitian_function_test(s[
    'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
fl.next('hermitian test')
fl.plot(residual)
print("best shift is",best_shift)
# {{{ slice out the FID appropriately and phase correct
# it
s.ft('t2')
s *= np.exp(1j*2*pi*best_shift*s.fromaxis('t2'))
s.ift('t2')
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
fl.next('signal vs. vd')
s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']].convolve('t2',10)#*-1 # bc inverted at high powers
s_forerror = s['ph2',coh_err['ph2']]['ph1',coh_err['ph1']]
# variance along t2 gives error for the mean, then average it across all other dimension, then sqrt for stdev
s_forerror.run(lambda x: abs(x)**2).mean_all_but(['vd']).run(sqrt)
s_sliced.mean('t2')#.set_error(s_forerror.data)
fl.plot(s_sliced,'o')
fl.plot(s_sliced.imag,'o')
fl.next('Spectrum - freq domain')
s = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']].convolve('t2',10)
fl.plot(s)
#s_sliced *= -1
fitfunc = lambda p, x: p[0]*(1-2*np.exp(-x*p[1]))
x = s_sliced.fromaxis('vd')
errfunc = lambda p: fitfunc(p,x).data - s_sliced.data.real
p_ini = [1.0,10.0]
p_opt,success = leastsq(errfunc, p_ini[:])
assert success in [1,2,3], "Fit did not succeed"
T1 = 1./p_opt[1]
print("T1:",T1,"s")
#new_x = nddata(r_[0:s_sliced.getaxis('vd')[-1]:500j],'vd')#.set_units('vd','s')
#fl.plot(fitfunc(p_opt,new_x), label='fit')
#gridandtick(gca())
s_sliced.name('data')
s_sliced.hdf5_write('proc_'+date+'_'+id_string+'_.h5')
fl.next('saved data')
fl.plot(s_sliced, 'o', label='data real')
fl.plot(s_sliced.imag, 'o', label='data imag')
fl.show();quit()

