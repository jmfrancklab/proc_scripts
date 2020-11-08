from pylab import *
from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
class flv(figlist_var):
    def show_complex(fl,d,**kwargs):
        fl.plot(abs(d),'k',alpha=0.5,**kwargs)
        fl.plot(d.real,'b',alpha=0.5,**kwargs)
        fl.plot(d.imag,'g',alpha=0.5,**kwargs)
        return
fl = flv()
for date,id_string,label_string,Rmax,whichslice in [
        ('200113','echo_DNP_TEMPOL_1','microwaves',1e2,('power',8)),
        ('191031','echo_5_mw_30dBm','+30 dBm microwaves',4e2,None),
        #('191031','echo_5_mw_34dBm','+34 dBm microwaves',4e2,None),
        #('191031','echo_5_mw_36dBm_2','+36 dBm microwaves',4e2,None),
        ]:
    fl.basenaem = id_string
    title_string = 'unenhanced'
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'ODNP_NMR_comp/old'))
    print(ndshape(s))
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1]).set_units('t2','s')
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.reorder('t2',first=False)
    if whichslice is not None:
        s = s[whichslice]
    s.ft(['ph1','ph2'])
    fl.next(id_string+'raw data - chunking coh')
    fl.image(s)
    s = s['ph1',1]['ph2',0]
    # {{{ slice only f-domain with signal
    s.ft('t2', shift=True)
    fl.next('frequency domain')
    fl.plot(abs(s), human_units=False)# no human units to make sure
    #                                   axvline matches
    slice_f = (-2e3,2e3)
    axvline(x=slice_f[0],
            color='k')
    axvline(x=slice_f[1],
            color='k')
    text(slice_f[1], abs(s).data.max()*0.75,
            ' (slice this region)',
            ha='left')
    s = s['t2':slice_f]
    s.ift('t2')
    # }}}
    fl.next('frequency filtered')
    fl.show_complex(s,human_units=False)
    pairs = s.contiguous(lambda x: abs(x) > abs(x.data).max()*0.5)
    longest_pair = diff(pairs).argmax()
    peak_location = pairs[longest_pair,:]
    axvline(x=peak_location[0],
            color='k')
    axvline(x=peak_location[1],
            color='k')
    text(peak_location[1], 0.85,
            ' (lines delineate peak of echo)',
            ha='left',
            transform=gca().transAxes)
    s.setaxis('t2',lambda x: x-peak_location.mean())
    # at this point, I need to get an axis that's in register with 0
    # -- this is what gives rise to "zero not in middle" below
    # I think I already have a function to do this
    s.register_axis({'t2':0})
    fl.next('crude centering')
    fl.show_complex(s)
    max_shift = diff(peak_location).item()/2
    dt = diff(s.getaxis('t2')[r_[0,1]]).item()
    max_shift = round(max_shift/dt)*dt
    fl.next('spectrum with crude centering')
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    fl.show_complex(s_sliced)
    s.ft('t2')
    # {{{ just check signs of shift
    s_left = s * exp(+1j*2*pi*max_shift*s.fromaxis('t2'))
    s_right = s * exp(-1j*2*pi*max_shift*s.fromaxis('t2'))
    s.ift('t2')
    s_left.ift('t2')
    s_right.ift('t2')
    fl.next('illustrate shifting')
    fl.plot(abs(s), label='original',alpha=0.5)
    fl.plot(abs(s_left), label='max left',alpha=0.5)
    fl.plot(abs(s_right), label='max right',alpha=0.5)
    # }}}
    shift_t = nddata(r_[-1:1:1000j]*max_shift, 'shift')
    t2_decay = exp(-s.fromaxis('t2')*nddata(r_[-Rmax:Rmax:1000j],'R2'))
    s_foropt = s.C
    s_foropt.ft('t2')
    s_foropt *= exp(1j*2*pi*shift_t*s_foropt.fromaxis('t2'))
    s_foropt.ift('t2')
    print("*** *** ***")
    print(ndshape(s_foropt))
    print("*** *** ***")
    s_foropt /= t2_decay
    s_foropt = s_foropt['t2':(-max_shift,max_shift)]
    print("max_shift",max_shift)
    dt = diff(s_foropt.getaxis('t2')[r_[0,1]]).item()
    # {{{ demand an odd number of points
    print(s_foropt.getaxis('t2')[r_[0,ndshape(s_foropt)['t2']//2,ndshape(s_foropt)['t2']//2+1,-1]])
    if ndshape(s_foropt)['t2'] % 2 == 0:
        print("chopped to get odd")
        s_foropt = s_foropt['t2',:-1]
    assert s_foropt.getaxis('t2')[s_foropt.getaxis('t2').size//2] == 0, 'zero not in the middle! -- does your original axis contain a 0?'
    # }}}
    # {{{ the middle point, at t=0 must be entirely real
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
    # replace following with time shift function
    s.ft('t2')
    s *= exp(1j*2*pi*best_shift*s.fromaxis('t2'))
    s.ift('t2')
    ph0 = s['t2':0.0]
    ph0 /= abs(ph0)
    s /= ph0
    fl.next('spectrum with optimized centering')
    s_sliced = s['t2':(0,None)].C
    s_sliced['t2',0] *= 0.5
    s_sliced.ft('t2')
    fl.show_complex(s_sliced)
    fl.next('Compare spectra')
    fl.plot(s_sliced.real,label='%s'%label_string)
    checking_time_slices = True
    print('best_R2',best_R2)
    for test_R2 in r_[best_R2]:# I have this as a list so we can mess around and choose a list of R2
        s_forplot = s / exp(-test_R2*s.fromaxis('t2'))
        fl.next(f'show superimposition R2={test_R2}', legend=True)
        fl.plot(abs(s_forplot), linewidth=3, alpha=0.5, label='abs')
        fl.plot(s_forplot.real, alpha=0.5, label='real')
        fl.plot(s_forplot.imag, alpha=0.5, label='imag')
        s_flipped = s_forplot.C.run(conj)
        s_flipped.setaxis('t2',lambda x: -x)
        fl.plot(abs(s_flipped), linewidth=3, alpha=0.5, label='abs flip')
        fl.plot(s_flipped.real, alpha=0.5, label='real flip')
        fl.plot(s_flipped.imag, alpha=0.5, label='imag flip')
fl.show()
