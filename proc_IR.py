from pyspecdata import *
from scipy.optimize import leastsq,minimize
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
# {{{ input parameters
date = '200203'
id_string = 'IR_pR_1_3'
clock_correction = 1.0829/998.253
nodename = 'signal'
filter_bandwidth = 1e3
# }}}
filename = date+'_'+id_string+'.h5'
s = nddata_hdf5(filename+'/'+nodename,
        directory = getDATADIR(exp_type = 'test_equip' ))
# it doesn't make sense to do any type of FT until the
# data is in the correct shape --- another reason why
# chunking should be done before the file is saved
print(s.getaxis('vd'))
s.chunk('t',['ph2','ph1','t2'],[4,2,-1]).set_units('t2','s')
s.setaxis('ph1',r_[0.:4.:2.]/4)
s.setaxis('ph2',r_[0.:4.]/4)
s.reorder(['ph1','ph2'])
fl.next('raw data')
s.ft('t2', shift=True)
s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
fl.image(s,interpolation='bicubic')
fl.next('filtered data')
s.ift('t2')
rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
s.setaxis(t2-rough_center)
s = s['t2':(-25e-3,25e-3)] # select only 50 ms in time domain, because it's so noisy
fl.image(s,interpolation='bicubic')
fl.next('before clock correction')
s.ft('t2')
fl.image(s)
fl.next('after clock correction')
s *= exp(-1j*s.fromaxis('vd')*clock_correction)
fl.image(s)
## My attempts to align the frequencies here do not work -- we should rather
## use the method from the YQ PNAS paper?
#fl.next('align frequencies')
#center_freqs = abs(s).argmax('t2')
#s.ift('t2')
#s *= exp(-1j*2*pi*center_freqs*s.fromaxis('t2'))
#s.ft('t2')
#fl.image(s)
fl.next('image coherence')
s.ft(['ph2','ph1'])
fl.image(s)
s.ift('t2')
residual,best_shift, best_R2 = hermitian_function_test(s['ph2',1]['ph1',0])
print("best shift is",best_shift)
# {{{ slice out the FID appropriately and phase correct
# it
s.ft('t2')
s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
s.ift('t2')
ph0 = s['t2',0]['ph2',1]['ph1',0]
if len(ph0.dimlabels) > 0:
    assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
    ph0 = zeroth_order_ph(ph0)
    print('phasing dimension as one')
else:
    ph0 = ph0/abs(ph0)
s /= ph0
fl.next('check phasing -- real')
fl.plot(s['ph2',1]['ph1',0])
gridandtick(gca())
fl.next('check phasing -- imag')
fl.plot(s['ph2',1]['ph1',0].imag)
gridandtick(gca())
s = s['t2':(0,None)]
s['t2',0] *= 0.5
fl.next('phased and FID sliced')
fl.image(s)
fl.next('phased and FID sliced -- frequency domain')
s.ft('t2')
# }}}
fl.image(s)
fl.next('signal vs. vd')
s_sliced = s['ph2',1]['ph1',0]
s_forerror = s['ph2',r_[0,2,3]]['ph1',1]
# variance along t2 gives error for the mean, then average it across all other dimension, then sqrt for stdev
s_forerror.run(lambda x: abs(x)**2).mean_all_but(['vd']).run(sqrt)
s_sliced.mean('t2').set_error(s_forerror.data)
fl.plot(s_sliced,'o')
fl.plot(s_sliced.imag,'o')
#min_vd = s.getaxis('vd')[abs(s).sum('t2').argmin('vd',raw_index=True).item()]
#est_T1 = min_vd/log(2)
#for x in range(len(s.getaxis('vd'))):
#    if s.getaxis('vd')[x] < min_vd:
#        s['vd',x] *= -1
#print("Estimated T1 is:",est_T1,"s")
fl.next('Spectrum - freq domain')
fl.plot(s['ph2',1]['ph1',0])
#fl.plot(s.reorder('vd',first=False), alpha=0.5)#, label='%s'%label_str)
##fl.next('Spectrum - waterfall')
##s_sliced.waterfall()
##s_sliced.sum('t2')
##s_sliced /= max(abs(s_sliced.data))
#fl.next('sum along t2')
#fl.plot(s.real,'o-')
#fl.plot(s.imag,'o-')
## perform fit for T1
#fl.next('T1 fit')
#x_data = s.getaxis('vd')
#y_data = s.data.real
#y_data /= max(abs(y_data))
#fl.plot(x_data,y_data,'o-',human_units=False)
#fitfunc = lambda p,x: p[0]*(1-2*exp(-x_data*p[1]))
#errfunc = lambda p_arg,x_arg,y_arg: fitfunc(p_arg,x_arg) - y_arg
#p_ini = [1.0,10.0]
#p1,success = leastsq(errfunc, p_ini[:], args=(x_data,y_data))
#assert success == 1, "Fit did not succeed"
#T1 = 1./p1[1]
#print("T1:",T1,"s")
#fl.show();quit()
#x_fit = linspace(x_data.min(),x_data.max(),5000)
#fl.next('ext')
#fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T1 = %0.2f ms)'%(T1*1e3), human_units=False)



fl.show();quit()

