from pylab import *
from pyspecdata import *
from proc_scripts import postproc_dict
zero_fill = False
with figlist_var() as fl:
    for filename,postproc,fslice,tslice,plen,max_kHz in [
            ('201228_Ni_sol_probe_nutation_amp_2','spincore_nutation_v2',
                (-28e3,23e3),(-2e-3,-1.3e-3),147e-6,150)
            ]:
        fl.basename = filename
        print('analyzing', filename)
        d = find_file(filename,exp_type='nutation',expno='nutation',postproc=postproc,
                lookup=postproc_dict,fl=fl)
        #d.chunk('t',['ph2','ph1','t2'],[2,4,-1])
        #d.set_units('t2','s')
        #d.ft('t2',shift=True)
        #fl.next('time domain')
        #fl.image(d.C.smoosh(['ph2','ph1'],'transient').reorder('transient').setaxis('transient','#').run(abs),
        #        interpolation='bilinear')
        #fl.show();quit()
        #d.reorder(['ph1','ph2'])
        #d.setaxis('ph2',r_[0:2]/4).setaxis('ph1',r_[0:4]/4)
        #if 'p_90' in d.dimlabels:
        #    d.set_units('p_90','s')
        #d.ft(['ph1','ph2'])
        #d = d['t2':(-20e3,25e3)]
        #fl.next('frequency domain')
        #fl.image(d)
        #fl.show();quit()
        print(d.get_prop('acq_params'))
        quit()
        d.ift('t2')
        print("max at",abs(d['ph1',1]['ph2',-2]).mean_all_but('t2').argmax('t2').item())
        d.setaxis('t2',lambda x: x-abs(d['ph1',1]['ph2',-2]).mean_all_but('t2').argmax('t2').item())
        fl.next('time domain--cropped log before slice')
        fl.image(d.C.cropped_log())
        #fl.show();quit()
        if tslice is not None:
            d = d['t2':tslice]
        print("signal max: %g"%abs(d['ph1',1]['ph2',-2].data.max()))
        fl.next('frequency domain after filter and t-slice')
        d.ft('t2')
        d = d['t2':fslice]
        fl.image(d)
        fl.show();quit()
        fl.next('slice out echo pathway')
        d = d['ph1',1]['ph2',-2]
        fl.image(d)
        #fl.show();quit()
        #gridandtick(gca(),gridcolor=[1,1,1])
        fl.next('FT')
        if zero_fill:
            d.ift('t2').ft('t2',pad=2**11)
        title('FT to get $\gamma B_1/a$')
        if 'amp' in d.dimlabels:
            d.setaxis('amp',lambda x: x*plen)
            d.set_units('amp','s')
            ind_dim = '\\tau_p a'
            d.rename('amp',ind_dim)
        if'p_90' in d.dimlabels:
            ind_dim = 'p_90'
        if zero_fill:
            d.ft(ind_dim, shift=True,pad=2**11)
        else:
            d.ft(ind_dim,shift=True)
        fl.image(d[ind_dim:(-1e3*max_kHz,1e3*max_kHz)])
        #fl.show();quit()
        fl.next('absFT')
        title('FT to get $\gamma B_1/a$')
        fl.image(abs(d[ind_dim:(-1e3*max_kHz,1e3*max_kHz)]))
        #gridandtick(gca(), gridcolor=[1,1,1])
        fl.show();quit()
                
