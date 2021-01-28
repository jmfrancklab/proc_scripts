from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import hermitian_function_test, zeroth_order_ph, recovery, integrate_limits, correl_align
from sympy import symbols, latex, Symbol
from matplotlib import *
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
coh_sel = {'ph1':0,
        'ph2':1}
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
fl.next('raw data -- coherence channels')
s.ft(['ph2','ph1'])
fl.image(s.C.setaxis('vd','#'))
s.ift('t2')
fl.next('time domain cropped log')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #').cropped_log())
#}}}
#{{{centering, hermitian function test and zeroth order phasing
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
    logger.info(strm('phasing dimension as one'))
else:
    logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
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
#}}}
#{{{slicing FID
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
#}}}
#{{{ Atttempting correlation alignment
s = s['t2':(-2.5e3,2.5e3)]
fl.next('before align')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
s.ift(['ph2','ph1'])
phasing = ndshape([4,2],['ph2','ph1']).alloc()
phasing.setaxis('ph1',r_[1,2]/4).setaxis('ph2',r_[0:4]/4)
phasing.ft(['ph2','ph1'])
phasing['ph2',1]['ph1',0] = 1
phasing.ift(['ph1','ph2'])
s /= phasing
fl.next('after FT coeff phase correction')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'), interpolation='bilinear')

s.smoosh(['vd','ph2','ph1'])
s.setaxis('vd','#')
s.reorder('t2',first=False)
fl.next('after smooshing in order')
fl.image(s)
fl.basename='first pass'
s = correl_align(s,align_phases=True,indirect_dim='vd',fl=fl)
fl.basename=None
s.ift('t2')
s.chunk('vd',['vd','ph1','ph2'],[-1,2,4])
s.setaxis('ph1',r_[0.,2.]/4)
s.setaxis('ph2',r_[0.,1.,2.,3.]/4)
s.setaxis('vd',r_[0:len(s.getaxis('vd'))])
s.reorder('ph1',first=True)
s.reorder('ph2',first=True)
s.reorder('t2',first=False)
s *= phasing
s.ft(['ph1','ph2'])
fl.next('time domain-after corr')
fl.image(s)
s.ft('t2')
fl.next('after alignment')
fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
fl.show();quit()
#}}}

fl.next('recovery curve')
s_signal = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
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
s = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
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
