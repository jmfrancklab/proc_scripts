from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
date = '200106'
label_str = 'none'
id_string0 = 'echo_off_1_6'
id_string1 = 'echo_on4_1'
filename0 = date+'_'+id_string0+'.h5'
filename1 = date+'_'+id_string1+'.h5'
nodename = 'signal'
s0 = nddata_hdf5(filename0+'/'+nodename,
        directory = getDATADIR(
            exp_type = 'test_equip'))
s1 = nddata_hdf5(filename1+'/'+nodename,
        directory = getDATADIR(
            exp_type = 'test_equip'))
nPoints = s0.get_prop('acq_params')['nPoints']
nEchoes = s0.get_prop('acq_params')['nEchoes']
nPhaseSteps = s0.get_prop('acq_params')['nPhaseSteps']
SW_kHz = s0.get_prop('acq_params')['SW_kHz']
nScans = s0.get_prop('acq_params')['nScans']
s0.reorder('t',first=True)
t2_axis = s0.getaxis('t')[0:nPoints/nPhaseSteps]
s0.setaxis('t',None)
s0.chunk('t',['ph2','ph1','t2'],[2,4,-1])
s0.setaxis('ph2',r_[0.,2.]/4)
s0.setaxis('ph1',r_[0.,1.,2.,3.]/4)
s0.setaxis('t2',t2_axis)
s0.setaxis('nScans',r_[0:nScans])
s0.reorder('t2',first=False)
s0.ft('t2',shift=True)
fl.next('raw data, chunked')
fl.image(abs(s0))
s0.ft(['ph1','ph2'])
fl.next('coherence')
fl.image(abs(s0))
s0 = s0['ph1',1]['ph2',0].C
s0.mean('nScans',return_error=False)
slice_f = (-3e3,3e3)
s0 = s0['t2':slice_f].C
s0.ift('t2')
max_data = abs(s0.data).max()
pairs = s0.contiguous(lambda x: abs(x) > max_data*0.5)
longest_pair = diff(pairs).argmax()
peak_location = pairs[longest_pair,:]
s0.setaxis('t2',lambda x: x-peak_location.mean())
s0.register_axis({'t2':0})
max_shift = diff(peak_location).item()/2
s0_sliced = s0['t2':(0,None)].C
s0_sliced['t2',0] *= 0.5
s0_sliced.ft('t2')
s0_ft = s0_sliced.C
fl.next('sliced')
fl.plot(s0_ft)
shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
t2_decay = exp(-s0.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
s0_foropt = s0.C
s0_foropt.ft('t2')
s0_foropt *= exp(1j*2*pi*shift_t*s0_foropt.fromaxis('t2'))
s0_foropt.ift('t2')
s0_foropt /= t2_decay
s0_foropt = s0_foropt['t2':(-max_shift,max_shift)]
print(s0_foropt.getaxis('t2'))
print(s0_foropt.getaxis('t2')[r_[0,ndshape(s0_foropt)['t2']//2,ndshape(s0_foropt)['t2']//2+1,-1]])
if ndshape(s0_foropt)['t2'] % 2 == 0:
    s0_foropt = s0_foropt['t2',:-1]
#assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2+1] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
ph0 = s0_foropt['t2':0.0]
ph0 /= abs(ph0)
s0_foropt /= ph0
s0_foropt /= max(abs(s0_foropt.getaxis('t2')))
# }}}
residual = abs(s0_foropt - s0_foropt['t2',::-1].runcopy(conj)).sum('t2')
residual.reorder('shift')
print(ndshape(residual))
minpoint = residual.argmin()
best_shift = minpoint['shift']
best_R2 = minpoint['R2']
s0.ft('t2')
s0 *= exp(1j*2*pi*best_shift*s0.fromaxis('t2'))
s0.ift('t2')
ph0 = s0['t2':0.0]
ph0 /= abs(ph0)
s0 /= ph0
s0_sliced = s0['t2':(0,None)].C
s0_sliced['t2',0] *= 0.5
s0_sliced.ft('t2')
fl.next('s0 Spectrum FT')
fl.plot(s0_sliced.real, alpha=0.5, label='%s'%label_str)


nPoints = s1.get_prop('acq_params')['nPoints']
nEchoes = s1.get_prop('acq_params')['nEchoes']
nPhaseSteps = s1.get_prop('acq_params')['nPhaseSteps']
SW_kHz = s1.get_prop('acq_params')['SW_kHz']
nScans = s1.get_prop('acq_params')['nScans']
s1.reorder('t',first=True)
t2_axis = s1.getaxis('t')[0:nPoints/nPhaseSteps]
s1.setaxis('t',None)
s1.chunk('t',['ph2','ph1','t2'],[2,4,-1])
s1.setaxis('ph2',r_[0.,2.]/4)
s1.setaxis('ph1',r_[0.,1.,2.,3.]/4)
s1.setaxis('t2',t2_axis)
s1.setaxis('nScans',r_[0:nScans])
s1.reorder('t2',first=False)
s1.ft('t2',shift=True)
fl.next('raw data, chunked')
fl.image(abs(s1))
s1.ft(['ph1','ph2'])
fl.next('coherence')
fl.image(abs(s1))
s1 = s1['ph1',1]['ph2',0].C
s1.mean('nScans',return_error=False)
slice_f = (-3e3,3e3)
s1 = s1['t2':slice_f].C
s1.ift('t2')
max_data = abs(s1.data).max()
pairs = s1.contiguous(lambda x: abs(x) > max_data*0.5)
longest_pair = diff(pairs).argmax()
peak_location = pairs[longest_pair,:]
s1.setaxis('t2',lambda x: x-peak_location.mean())
s1.register_axis({'t2':0})
max_shift = diff(peak_location).item()/2
s1_sliced = s1['t2':(0,None)].C
s1_sliced['t2',0] *= 0.5
s1_sliced.ft('t2')
s1_ft = s1_sliced.C
fl.next('sliced')
fl.plot(s1_ft)
shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
t2_decay = exp(-s1.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
s1_foropt = s1.C
s1_foropt.ft('t2')
s1_foropt *= exp(1j*2*pi*shift_t*s1_foropt.fromaxis('t2'))
s1_foropt.ift('t2')
s1_foropt /= t2_decay
s1_foropt = s1_foropt['t2':(-max_shift,max_shift)]
print(s1_foropt.getaxis('t2'))
print(s1_foropt.getaxis('t2')[r_[0,ndshape(s1_foropt)['t2']//2,ndshape(s1_foropt)['t2']//2+1,-1]])
if ndshape(s1_foropt)['t2'] % 2 == 0:
    s1_foropt = s1_foropt['t2',:-1]
#assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2+1] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
ph0 = s1_foropt['t2':0.0]
ph0 /= abs(ph0)
s1_foropt /= ph0
s1_foropt /= max(abs(s1_foropt.getaxis('t2')))
# }}}
residual = abs(s1_foropt - s1_foropt['t2',::-1].runcopy(conj)).sum('t2')
residual.reorder('shift')
print(ndshape(residual))
minpoint = residual.argmin()
best_shift = minpoint['shift']
best_R2 = minpoint['R2']
s1.ft('t2')
s1 *= exp(1j*2*pi*best_shift*s1.fromaxis('t2'))
s1.ift('t2')
ph0 = s1['t2':0.0]
ph0 /= abs(ph0)
s1 /= ph0
g_t = s1/s0
#g_t = g_t['t2':(0,0.004)]
fl.next('idk')
fl.plot(g_t)
g_t.ft('t2')
fl.next('f')
fl.plot(g_t)

fl.show();quit()
