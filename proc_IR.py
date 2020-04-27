from pyspecdata import *
from scipy.optimize import leastsq,minimize
from pyspecdata.load_files.bruker_nmr import bruker_data
from proc_scripts import hermitian_function_test, zeroth_order_ph,load_data
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
# {{{ input parameters
#clock_correction = 1.785
#clock_correction = 0
nodename = 'signal'
filter_bandwidth = 5e3
coh_sel = {'ph1':0,
        'ph2':1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
expno = 2
s = load_data('w12_200309')
# }}}
#for exp_name,expno in [
#        ('w12_200309',2),
#        ]:
#    s = find_file(exp_name, exp_type='NMR_Data_AG', dimname = 'indirect', expno=expno)
s.chunk('indirect',['indirect','ph1','ph2'],[-1,4,2])
s.reorder(['ph2','ph1']).set_units('t2','s')
fl.next('raw data')
fl.image(s)
fl.next('after clock correction')
s.rename('indirect','vd')
#s *= exp(-1j*s.fromaxis('vd')*clock_correction)
#fl.image(s)
fl.next('raw data -- coherence channels')
s.ft(['ph2','ph1'])
fl.image(s)
fl.next('filtered + rough centered data')
s.ft('t2', shift=True)
s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
s.ift('t2')
rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
s.setaxis(t2-rough_center)
#s = s['t2':(-25e-3,25e-3)] # select only 50 ms in time domain, because it's so noisy
s.ft('t2')
#{{{ aligning freq
## My attempts to align the frequencies here do not work -- we should rather
## use the method from the YQ PNAS paper?
#fl.next('align frequencies')
#center_freqs = abs(s).argmax('t2')
#s.ift('t2')
#s *= exp(-1j*2*pi*center_freqs*s.fromaxis('t2'))
#s.ft('t2')
#fl.image(s)
#}}}
s.ift('t2')
residual,best_shift = hermitian_function_test(s[
    'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
fl.next('hermitian test')
fl.plot(residual)
print("best shift is",best_shift)
# {{{ slice out the FID appropriately and phase correct
# it
s.ft('t2')
s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
s.ift('t2')
fl.next('time domain after hermitian test')
fl.image(s)
ph0 = s['t2':0]['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
print(ndshape(ph0))
if len(ph0.dimlabels) > 0:
    assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
    ph0 = zeroth_order_ph(ph0, fl=fl)
    print('phasing dimension as one')
else:
    print("there is only one dimension left -- standard 1D zeroth order phasing")
    ph0 = ph0/abs(ph0)
s /= ph0
fl.next('frequency domain -- after hermitian function test and phasing')
s.ft('t2')
fl.image(s)

s.ift('t2')
fl.next('check phasing -- real')
fl.plot(s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']])
#gridandtick(gca())
fl.next('check phasing -- imag')
fl.plot(s[
    'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']].imag)

#gridandtick(gca())
s = s['t2':(0,None)]
s['t2',0] *= 0.5
fl.next('phased and FID sliced')
fl.image(s)
fl.next('phased and FID sliced -- frequency domain')
s.ft('t2')
# }}}
fl.image(s)
fl.next('signal vs. vd')
s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]*-1 # bc inverted at high powers
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
print("T1:",T1,"s")
fl.next('fit')
fl.plot(s_sliced, 'o', label='data')
fl.plot(s_sliced.imag, 'o', label='data')
new_x = nddata(r_[0:s_sliced.getaxis('vd')[-1]:500j],'vd')#.set_units('vd','s')
fl.plot(fitfunc(p_opt,new_x), label='fit')
gridandtick(gca())
fl.show();quit()
