from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import hermitian_function_test, zeroth_order_ph, recovery, integrate_limits, correl_align, ph1_real_Abs
from sympy import symbols, latex, Symbol
from matplotlib import *
from scipy.signal import tukey
import matplotlib.pyplot as plt
import numpy as np
from sympy import exp as s_exp
init_logging('debug')
fl = figlist_var()
t2 = symbols('t2')
# {{{ input parameters
date = '210127'
id_string = 'OHTEMPO10mM_cap_probe_IR_1'
nodename = 'signal'
filter_bandwidth = 5e3
ph1_val = 0
ph2_val = 1
coh_err = {'ph1':1,# coherence pathways to use for error -- note that this
        #             should ideally be pathways that do NOT include any known
        #             artifacts
        'ph2':r_[0,2,3]}
# }}}
#{{{preprocessing we can later add to postproc_dict
filename = date+'_'+id_string+'.h5'
s = nddata_hdf5(filename+'/'+nodename,
        directory = getDATADIR(exp_type='test_equip'))
vd_axis = s.getaxis('vd')
s.reorder(['ph2','ph1']).set_units('t2','s')
s.ft('t2',shift=True)
fl.next('raw data')
fl.image(s.C.setaxis('vd','#'))
fl.next('freq domain -- cropped log')
s.ft(['ph2','ph1'])
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #').cropped_log())
f_range = (-4000,4000)
s = s['t2':f_range]
s.ift('t2')
fl.next('time domain cropped log')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #').cropped_log())
t_range = (0,35e-3)
s = s['t2':t_range]
#}}}
#{{{centering, hermitian function test and zeroth order phasing
rough_center = abs(s).convolve('t2',10).mean_all_but('t2').argmax('t2').item()
s.setaxis(t2-rough_center)
fl.next('rough centering')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
apod = s.fromaxis('t2')* 0 
apod['t2':(None,500)] = apod['t2':(None,500)].run(lambda x:tukey(len(x)))
s *= apod
#}}}
#{{{phasing the aligned data
#s.ift('t2')
best_shift = hermitian_function_test(s[
    'ph2',1]['ph1',0].C.mean('vd'))
logger.info(strm("best shift is", best_shift))
s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
fl.next('time domain after hermitian test')
fl.image(s.C.setaxis(
'vd','#').set_units('vd','scan #'))

ph0 = s['t2':0]['ph2',ph2_val]['ph1',ph1_val]
if len(ph0.dimlabels) > 0:
    assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
    ph0 = zeroth_order_ph(ph0)
    logger.info(strm('phasing dimension as one'))
else:
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
    ph0 = ph0/abs(ph0)
s /= ph0
dw = np.diff(s.getaxis('t2')[r_[0,-1]]).item()
print("THIS IS SHAPE OF S",ndshape(s))
s = ph1_real_Abs(s,dw,fl=fl)
fl.show();quit()
fl.next('phased data time domain')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
fl.next('phased data freq domain')
s.ft('t2')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
fl.show();quit()
#}}}
#{{{slicing FID
#s = s['t2':(0,None)]
#s['t2',0] *= 0.5
#fl.next('phased and FID sliced')
#fl.image(s.C.setaxis(
#'vd','#').set_units('vd','scan #'))
#fl.next('phased and FID sliced -- frequency domain')
#s.ft('t2')
#fl.image(s.C.setaxis(
#'vd','#').set_units('vd','scan #'))
#}}}
s.ift('t2')
#{{{ Atttempting correlation alignment
fl.next('before aligning prior to null')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
for j in range(ndshape(s)['vd']):
    if s['vd',j].C.sum('t2')<1:
        s['vd',j] *= -1
s_final = s.C * 0
s.ft('t2')
fl.next('before splitting and correl align')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
s_before = s['vd',:2]
s_after = s['vd',3:]
fl.next('prior to alignment s_before- coherence dimension')
fl.image(s_before)
fl.next('prior to alignment s_after- coherence dimension')
fl.image(s_after.C.setaxis('vd','#').set_units('vd','scan #'))
s_before.ift(['ph2','ph1'])
s_before.ift('t2')
s_before.ft('t2')
fl.next('prior to alignment')
fl.image(s_before, interpolation = 'bilinear')
fl.basename='first pass'
#opt_shift_after,sigma_after = correl_align(s_after,indirect_dim='vd',ph1_selection=ph1_val,ph2_selection=ph2_val,fl=fl)
opt_shift_before,sigma_before = correl_align(s_before,indirect_dim='vd',ph1_selection=ph1_val,ph2_selection=ph2_val,fl=fl)
fl.basename= None
s_before.ift('t2')
s_before *= np.exp(-1j*2*pi*opt_shift_before*s_before.fromaxis('t2'))
#s_after *= np.exp(-1j*2*pi*opt_shift_after*s_after.fromaxis('t2'))
fl.next('out of correl alignment')
fl.image(s_before)
s_before.reorder('ph1',first=True)
s_before.reorder('ph2',first=True)
s_before.reorder('t2',first=False)
s_before.ft(['ph1','ph2'])

s_after.ift(['ph2','ph1'])
s_after.ift('t2')
s_after.ft('t2')
fl.next('prior to alignment')
fl.image(s_after.C.setaxis('vd','#').set_units('vd','scan #'), interpolation = 'bilinear')
fl.basename='first pass'
opt_shift_after,sigma_after = correl_align(s_after,indirect_dim='vd',ph1_selection=ph1_val,ph2_selection=ph2_val,fl=fl)
fl.basename= None
s_after.ift('t2')
s_after *= np.exp(-1j*2*pi*opt_shift_after*s_after.fromaxis('t2'))
fl.next('out of correl alignment')
fl.image(s_after.C.setaxis('vd','#').set_units('vd','scan #'))
s_after.reorder('ph1',first=True)
s_after.reorder('ph2',first=True)
s_after.reorder('t2',first=False)
s_after.ft(['ph1','ph2'])
s_final['vd',:2] = s_before
s_final['vd',3:] = s_after
s = s_final

fl.next('coherence domain-after corr')
s.ft('t2')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
s.ift('t2')
s *= -1
fl.next('after alignment')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
#}}}

fl.next('recovery curve')
s_signal = s['ph2',ph2_val]['ph1',ph1_val]
# {{{ here we use the inactive coherence pathways to determine the error
#     associated with the data
s_forerror = s['ph2',coh_err['ph2']]['ph1',coh_err['ph1']]
# variance along t2 gives error for the mean, then average it across all other dimension, then sqrt for stdev
logger.info(strm(ndshape(s_forerror)))
frq_slice = integrate_limits(s_signal)
s_signal = s_signal['t2':frq_slice]
s_forerror = s_forerror['t2':frq_slice]
s_forerror.run(lambda x:abs(x)**2).mean_all_but(['vd','t2']).integrate('t2').run(sqrt)
s_signal.integrate('t2').set_error(s_forerror.data)
# }}}
logger.info(strm("here is what the error looks like",s_signal.get_error()))
fl.plot(s_signal,'o',label='real')
fl.plot(s_signal.imag,'o',label='imaginary')
fl.next('Spectrum - freq domain')
s = s['ph2',ph2_val]['ph1',ph1_val]
fl.plot(s)
fitfunc = lambda p, x: p[0]*(1-2*np.exp(-x*p[1]))
x = s_signal.fromaxis('vd')
errfunc = lambda p: fitfunc(p,x).data - s_signal.data.real
p_ini = [1.0,10.0]
#p_opt,success = leastsq(errfunc, p_ini[:])
#assert success in [1,2,3], "Fit did not succeed"
#T1 = 1./p_opt[1]
#logger.info(strm("T1:",T1,"s"))
f = fitdata(s_signal)
error = fitdata(s_forerror)
M0,Mi,R1,vd = symbols("M_0 M_inf R_1 vd", real=True)
error.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
f.functional_form = Mi + (M0-Mi)*s_exp(-vd*R1)
error.fit()
f.fit()
logger.info(strm("output:",f.output()))
logger.info(strm("latex:",f.latex()))
T1 = 1./f.output('R_1')
fl.next('fit')
fl.plot(s_signal,'o', label='actual data')
fl.plot(s_signal.imag,'o',label='actual imaginary')
fl.plot(f.eval(100),label='fit')
fl.plot(error.eval(100),label='fit error')
plt.text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes,size='large',
        horizontalalignment='center',color='k')
ax = plt.gca()
fl.show();quit()

