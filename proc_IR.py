from pyspecdata import *
from scipy.optimize import leastsq,minimize
fl = figlist_var()
date = '200203'
id_string = 'IR_pR_1_3'
filename = date+'_'+id_string+'.h5'
nodename = 'signal'
s = nddata_hdf5(filename+'/'+nodename,
        directory = getDATADIR(exp_type = 'test_equip' ))
s.rename('t','t2').set_units('t2','s')
s.ft('t2',shift=True)
clock_correction = 1.0829/998.253
s *= exp(-1j*s.fromaxis('vd')*clock_correction)
fl.next('raw data - clock correction')
fl.image(s['t2':(-1e3,1e3)])
s.ift('t2')
s.chunk('t2',['ph2','ph1','t2'],[4,2,-1])
#s.setaxis('t2',t2_axis[nPoints])
s.setaxis('ph1',r_[0.,2.]/4)
s.setaxis('ph2',r_[0.,1.,2.,3.]/4)
fl.next('image')
s.reorder('vd',first=False)
s.reorder('t2',first=False)
s.ft(['ph2','ph1'])
fl.image(s['t2':(-1e3,1e3)].C.ft('t2'))
s = s['ph2',1]['ph1',0]
slice_f = (-600,800)
s.ft('t2')
s = s['t2':slice_f]
s.ift('t2')
fl.next('rough center')
max_data = abs(s.data).max()
vd_pairs = []
for x in xrange(len(s.getaxis('vd'))):
    pairs = s['vd',x].contiguous(lambda t: abs(t) > max_data*0.5)
    try:
        longest_pair = diff(pairs).argmax()
        peak_location = pairs[longest_pair,:]
        print peak_location
        # deciding to find phase corrections from first peak
        #if x == len(s.getaxis('vd'))-1:
        if x == 0:
            max_shift = diff(peak_location).item()/2
        vd_pairs.append(peak_location.mean())
    except:
        print "No max found"
        pass
peak_location = array(vd_pairs).mean()
s.setaxis('t2',lambda x: x-peak_location)
#s.register_axis({'t2':0})
fl.image(s.C.ft('t2'))
s_sliced = s['t2':(0,None)]
s_sliced['t2',0] *= 0.5
fl.next('all FIDs')
fl.plot(s_sliced.C.reorder('vd',first=False))
fl.next('all FIDs - FT')
fl.plot(s_sliced.C.ft('t2').reorder('vd',first=False))
shift_t = nddata(r_[-1:1:200j]*max_shift, 'shift')
#s_foropt = s['vd',-1].C
s_foropt = s['vd',0].C
s_foropt.ft('t2')
s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
s_foropt.ift('t2')
s_foropt = s_foropt['t2':(-max_shift,max_shift)]
if ndshape(s_foropt)['t2'] % 2 == 0:
    s_foropt = s_foropt['t2',:-1]
ph0 = s_foropt['t2':0.0]
ph0 /= abs(ph0)
s_foropt /= ph0
s_foropt /= max(abs(s_foropt.getaxis('t2')))
residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
residual.reorder('shift')
minpoint = residual.argmin()
best_shift = minpoint['shift']
s.ft('t2')
s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
s.ift('t2')
ph0 = s['t2':0.0]
ph0 /= abs(ph0)
s /= ph0
fl.next('phased')
fl.image(s.C.ft('t2'))
s_sliced = s['t2':(0,None)]
s_sliced['t2',0] *= 0.5
s_sliced.ft('t2')
min_vd = s_sliced.getaxis('vd')[abs(s_sliced).sum('t2').argmin('vd',raw_index=True).item()]
est_T1 = min_vd/log(2)
for x in range(len(s_sliced.getaxis('vd'))):
    if s_sliced.getaxis('vd')[x] < min_vd:
        s_sliced['vd',x] *= -1
s_sliced.sum('t2')
s_sliced /= max((s_sliced.data))
fl.next('sum along t2 - real')
fl.plot(s_sliced.real,'o-')
fl.plot(s_sliced.imag,'o-')
fl.show();quit()
