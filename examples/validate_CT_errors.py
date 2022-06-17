from pyspecdata import *
from pylab import *
from matplotlib import *
from pyspecProcScripts import *
from pyspecProcScripts.correlation_alignment import correl_align
import numpy as np
import sympy as s
from collections import OrderedDict
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery
from numpy.random import seed
seed(120)
# sphinx_gallery_thumbnail_number = 4

fl = figlist_var()
signal_pathway = {"ph1": 1, "ph2": 0}
excluded_pathways = [(0, 0), (0, 3)]
colors = ["r", "darkorange", "gold", "g", "c", "b", "m", "lightcoral"]
# {{{ generate the fake data
t2, nScans, ph1, ph2 = s.symbols('t2 nScans ph1 ph2')
data = fake_data(
    (23*(1+1e-8*nScans)*(s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi))),
    OrderedDict([
        ("nScans" , nddata(r_[0.0:10000], "nScans")),
        ("ph2" , nddata(r_[0.0,2.0] / 4, "ph2")),
        ("ph1" , nddata(r_[0., 1., 2., 3.] / 4, "ph1")),
        ("t2" , nddata(r_[0:0.085:256j], "t2"))]),
        signal_pathway, scale = 15.0)
# {{{ just have the data phased
data.labels({'ph2':r_[0.0,2.0]/4,'ph1':r_[0.0,1.0,2.0,3.0]/4})
data.reorder(["ph1", "ph2"])
data.ft('t2')
data /= sqrt(ndshape(data)['t2'])*data.get_ft_prop('t2','dt')
data.ift('t2')
# }}}
data = data["t2" : (0,None)]
data["t2":0] *= 0.5
data.ft("t2")
fl.next('Fake Data -- 0.0 scale for freq var')
fl.image(data)
# }}}
# {{{Normalization
int_frq_slice_start = integrate_limits(select_pathway(data['nScans',0].C, signal_pathway),cutoff = 0.1)
int_frq_slice_end = integrate_limits(select_pathway(data['nScans',-1].C,signal_pathway),cutoff = 0.1)
int_slice = (int_frq_slice_start[0],int_frq_slice_end[-1])
spacing = int_slice[-1]-int_slice[0]
offres_start = int_slice[-1]+500
s_integral = select_pathway(data['t2':int_slice].C, signal_pathway).integrate('t2')
avg_d = s_integral.C.mean().real.item()
s_integral /= avg_d
data /= avg_d
# }}}
integral_diff_sq=abs(s_integral.C.real - s_integral.C.mean().real)**2
f= data.C.getaxis('t2')
df = f[1]-f[0]
fl.next('limits')
fl.plot(select_pathway(data.C.mean('nScans'),signal_pathway))
plt.axvline(int_slice[0],color='blue')
plt.axvline(int_slice[-1],color='blue')
#{{{control var
s_intregion = s_integral.real
s_intregion.run(var,'nScans')
#}}}
#{{{ inactive pathways
s_inactive = data['t2':int_slice].C
s_inactive['ph1',1]['ph2',0] = 0
N = ndshape(s_inactive.C)['t2']
s_inactive.run(var,'t2')
s_inactive /= 2 # equiv to variance of real datapoints
s_inactive *= df **2
s_inactive *= N
#}}}
#{{{off resonance slice
s_offres = select_pathway(data['t2':(offres_start,None)].C,signal_pathway)
s_offres = s_offres['t2',0:N]
assert (ndshape(s_offres)['t2']) == N
plt.axvline(s_offres.C.getaxis('t2')[0],color='red')
plt.axvline(s_offres.C.getaxis('t2')[-1],color='red')
s_offres = s_offres.real # only keep the real part
s_offres.run(var,'t2')
# {{{ convert from error of datapoints to propagated
#     error for integral
s_offres *= df **2
s_offres *= N
s_offres_avg = s_offres.C.mean().item()
#}}}
# }}}
fl.next('variances')
plt.axhline(y = s_intregion.data, color='k',label='variance of integral')
fl.plot(s_offres,'o',color='red',label = 'variance from off res slice')
plt.axhline(y = s_offres_avg,color='red', label = 'averaged off res')
s_inactive.smoosh(['ph1','ph2'])
for j in range(ndshape(s_inactive)['ph1']):
    fl.plot(s_inactive['ph1',j],'o', label='variance of ph1 = %d'%j)
avg_inactive = (1/(ndshape(s_inactive.C)['ph1'])*s_inactive.C.sum('ph1')).mean().item()
plt.axhline(y = avg_inactive,color = 'green',label='averaged inactive')
fl.plot(integral_diff_sq,'o',color='k',label = 'integral_diff_sq')
print("variance of integral",s_intregion)
print("variance of integral, predicted based on propagation of offres",s_offres_avg)
#print("variance of integral, predicted based on propagation of inactive:",s_inactive)
#for j in range(ndshape(s_inactive)['ph1']):
#    print(s_inactive.getaxis('ph1')[j],'-->',s_inactive['ph1',j])
# in the following the denominator doesn't include the
# active pathway
print('and averaged:', avg_inactive)
        #1/(ndshape(s_inactive)['ph1']-1)*s_inactive.sum('ph1'))
fl.show()
