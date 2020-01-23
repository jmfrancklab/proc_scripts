from pyspecdata import *
from scipy.optimize import leastsq,minimize
fl = figlist_var()
date = '200122'
id_string = 'IR_water_6'
filename = date+'_'+id_string+'.h5'
nodename = 'signal'
s = nddata_hdf5(filename+'/'+nodename,
        directory = getDATADIR(exp_type = 'test_equip' ))
s.rename('t','t2').set_units('t2','s')
print s.getaxis('vd')
fl.next('raw data - no clock correction')
fl.image(s)
s.ft('t2',shift=True)
clock_correction = 0 # radians per second
clock_correction = 1.0829/998.253
#clock_correction = -0.58434/9.88461 # for IR_1
#clock_correction = -4.33/9.72
s *= exp(-1j*s.fromaxis('vd')*clock_correction)
s.ift('t2')
fl.next('raw data - clock correction')
fl.image(s)
nPoints = 1024*2
t2_axis = s.getaxis('t2')
s.setaxis('t2',None)
s.chunk('t2',['ph2','ph1','t2'],[4,2,-1])
s.setaxis('t2',t2_axis[nPoints])
s.setaxis('ph1',r_[0.,2.]/4)
s.setaxis('ph2',r_[0.,1.,2.,3.]/4)
print ndshape(s)
fl.next('image')
fl.image(s)
s.ft(['ph2','ph1'])
s.ft('t2')
fl.next('image coherence')
fl.image(s['vd',0])
s = s['ph2',1]['ph1',0].C
slice_f = (-3e3,3e3)
s = s['t2':slice_f].C
s.ift('t2')
max_data = abs(s.data).max()
pairs = s['vd',0].contiguous(lambda x: abs(x) > max_data*0.5)
longest_pair = diff(pairs).argmax()
peak_location = pairs[longest_pair,:]
s.setaxis('t2',lambda x: x-peak_location.mean())
s.register_axis({'t2':0})
max_shift = diff(peak_location).item()/2
s_sliced = s['t2':(0,None)].C
s_sliced['t2',0] *= 0.5
s_ft = s_sliced.C
fl.next('all corrected FIDs')
fl.plot(s_ft.reorder('vd',first=False))
shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:200j],'R2'))
s_foropt = s.C
s_foropt.ft('t2')
s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
s_foropt.ift('t2')
s_foropt /= t2_decay
s_foropt = s_foropt['t2':(-max_shift,max_shift)]
print s_foropt.getaxis('t2')
print s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]]
if ndshape(s_foropt)['t2'] % 2 == 0:
    s_foropt = s_foropt['t2',:-1]
ph0 = s_foropt['t2':0.0]
ph0 /= abs(ph0)
s_foropt /= ph0
s_foropt /= max(abs(s_foropt.getaxis('t2')))
residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
residual.reorder('shift')
print ndshape(residual)
minpoint = residual.argmin()
best_shift = minpoint['shift']
best_R2 = minpoint['R2']
s.ft('t2')
s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
s.ift('t2')
ph0 = s['t2':0.0]
ph0 /= abs(ph0)
s /= ph0
s_sliced = s['t2':(0,None)].C
s_sliced['t2',0] *= 0.5
s_sliced.ft('t2')
min_vd = s_sliced.getaxis('vd')[abs(s_sliced).sum('t2').argmin('vd',raw_index=True).item()]
est_T1 = min_vd/log(2)
for x in xrange(len(s_sliced.getaxis('vd'))):
    if s_sliced.getaxis('vd')[x] < min_vd:
        s_sliced['vd',x] *= -1
print "Estimated T1 is:",est_T1,"s"
fl.next('Spectrum - freq domain')
fl.plot(s_sliced.reorder('vd',first=False), alpha=0.5)#, label='%s'%label_str)
fl.next('Spectrum - waterfall')
s_sliced.waterfall()
s_sliced.sum('t2')
fl.next('sum along t2')
fl.plot(s_sliced.real,'o-')
fl.plot(s_sliced.imag,'o-')
# perform fit for T1
fl.next('T1 fit')
x_data = s_sliced.getaxis('vd')
y_data = s_sliced.data.real
y_data /= max(abs(y_data))
fl.plot(x_data,y_data,'o-',human_units=False)
fitfunc = lambda p,x: p[0]*(1-2*exp(-x_data*p[1]))
errfunc = lambda p_arg,x_arg,y_arg: fitfunc(p_arg,x_arg) - y_arg
p_ini = [1.0,10.0]
p1,success = leastsq(errfunc, p_ini[:], args=(x_data,y_data))
assert success == 1, "Fit did not succeed"
T1 = 1./p1[1]
print "T1:",T1,"s"
fl.show();quit()
x_fit = linspace(x_data.min(),x_data.max(),5000)
fl.next('ext')
fl.plot(x_fit, fitfunc(p1, x_fit),':', label='fit (T1 = %0.2f ms)'%(T1*1e3), human_units=False)



fl.show();quit()

