from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
# {{{ in this section, I define some manual manipulations that might be
# specific to the particular type of data, with the idea being that these could
# be modified for different acquisition conditions / pulse sequences
# need pyspecdata from 11/22/19 or sooner
class nddata_withopt (nddata):
    def __init__(self,d):
        super(self.__class__,self).__init__(0)
        self.set_to(d)
        s.coh_procd = False
    def chunk_dims(s):
        "divide up the time dimension appropriately"
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1]).set_units('t2','s')
    def initial_filter(self):
        "FT but also throw out f-domain data that's just noise"
        self.ft('t2',shift=True)
        self.set_to(self['t2':(-1e3,1e3)])
        return self
    def proc_coherence(s):
        s.coh_procd = True
        s.setaxis('ph2',r_[0:4:2]/4.).setaxis('ph1',r_[0:4]/4.)
        s.ft(['ph2','ph1'])
    def phcyc(s):
        if not s.coh_procd:
            s.proc_coherence()
        s.set_to(s['ph1',1]['ph2',0])
    def return_timeordered(s):
        "to see if we have frequency drift, give the data in time order"
        temp = abs(s)
        temp.reorder('power') # first/outer loop
        temp.smoosh(['power','ph2','ph1'],dimname='indirect')
        temp.reorder('t2',first=False)
        temp.setaxis('indirect',None)
        return temp
# }}}
fl = figlist_var()
def get_W(dBm):
    return 10**(dBm/10.)*1e-3
dBm_list = [0., 30., 34., 36.]
W_list = ones_like(dBm_list)
for x in xrange(4):
    W_list[x] = get_W(dBm_list[x])
enhancement = []
find_phase_params = True # phase params found for first dataset will be applied
                         # to all subsequently processed datasets
datalist = []
mw_power_list = []
for date,id_string,label_string,mw_power in [
        ('191031','echo_5_4','no microwaves',-9999),
        ('191031','echo_5_mw_30dBm','+30 dBm microwaves',30),
        ('191031','echo_5_mw_34dBm','+34 dBm microwaves',34),
        ('191031','echo_5_mw_36dBm_2','+36 dBm microwaves',36),
        ]:
    title_string = 'unenhanced'
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    datalist.append(s)
    mw_power_list.append(mw_power)
# {{{ this is the main routine -- details of how this
#     works are given by the nddata_withopt class
s = concat(datalist,'power').setaxis('power',array(mw_power_list))
s = nddata_withopt(s)
s.chunk_dims()
s.reorder('t2',first=False)
fl.next('raw signal')
fl.image(s)
fl.next('FT of raw, filtered signal')
s.initial_filter()
fl.image(abs(s))
fl.next('data in time order')
fl.image(s.return_timeordered())
fl.next('look at coherence purity')
s.proc_coherence()
fl.image(s)
fl.next('abs time domain')
s.phcyc()
s.ift('t2', pad=2048) # just to make it pretty
fl.plot(s, alpha=0.5, label_format_string='%f dBm')
fl.plot(s.imag, alpha=0.2)
fl.show()
# }}}
if False:
    fl.next('frequency domain')
    fl.plot(abs(s))
    slice_f = (-1e3,1e3)
    s = s['t2':slice_f]
    s.ift('t2')
    max_data = abs(s.data).max()
    pairs = s.contiguous(lambda x: abs(x) > max_data*0.5)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    s.setaxis('t2',lambda x: x-peak_location.mean())
    s.register_axis({'t2':0})
    max_shift = diff(peak_location).item()/2
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    if find_phase_params:
        shift_t = nddata(r_[-1:1:1000j]*max_shift, 'shift')
        t2_decay = exp(-s.fromaxis('t2')*nddata(r_[0:1e3:1000j],'R2'))
        s_foropt = s.C
        s_foropt.ft('t2')
        s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
        s_foropt.ift('t2')
        s_foropt /= t2_decay
        s_foropt = s_foropt['t2':(-max_shift,max_shift)]
        print s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]]
        if ndshape(s_foropt)['t2'] % 2 == 0:
            s_foropt = s_foropt['t2',:-1]
        assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2+1] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
        ph0 = s_foropt['t2':0.0]
        ph0 /= abs(ph0)
        s_foropt /= ph0
        s_foropt /= max(abs(s_foropt.getaxis('t2')))
        # }}}
        residual = abs(s_foropt - s_foropt['t2',::-1].runcopy(conj)).sum('t2')
        fl.next('cost function')
        residual.reorder('shift')
        fl.image(residual)
        fl.plot(residual.C.argmin('shift').name('shift'),'x')
        minpoint = residual.argmin()
        best_shift = minpoint['shift']
        best_R2 = minpoint['R2']
        fl.plot(best_R2,best_shift,'o')
        find_phase_params = False
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    ph0 = s['t2':0.0]
    ph0 /= abs(ph0)
    s /= ph0
    fl.next('Aer spectra')
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    fl.plot(s_sliced.real, alpha=0.5, label='%s'%label_string)
    enhancement.append(s_sliced.real.sum('t2').item())
    fl.next('E(p)')
    fl.plot(W_list,array(enhancement),'o-')
    fl.show('aer_prelim_ODNP_2.pdf')
