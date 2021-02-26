from pylab import *
from pyspecdata import *
from proc_scripts import postproc_dict
zero_fill = False
with figlist_var() as fl:
    for filename,postproc,fslice,tslice,max_kHz in [
            ('210120_Ni_sol_probe_nutation_amp_1','spincore_nutation_v2',
                (-15e3,15e3),(-0.75e-3,0.75e-3),200)
            ]:
        
        fl.basename = filename
        print('analyzing', filename)
        d = find_file(filename,exp_type='nutation',expno='nutation',postproc=postproc,
                lookup=postproc_dict,fl=fl)
        #print(d.get_prop('acq_params'))
        plen = d.get_prop('acq_params')['p90_us']
        plen *= 10**-6
        d = d['t2':(-20e3,20e3)]
        fl.next('frequency domain--after slice')
        fl.image(d)
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
        #fl.show();quit()
        fl.next('slice out echo pathway')
        d = d['ph1',1]['ph2',-2]
        if 'amp' in d.dimlabels:
            d.setaxis('amp',lambda x: x*plen)
            d.set_units('amp','s')
            ind_dim = '\\tau_p a'
            d.rename('amp',ind_dim)
        if'p_90' in d.dimlabels:
            ind_dim = 'p_90'
        fl.image(d)
        fl.next('time domain')
        d.ift('t2')
        fl.image(d)
        d.ft('t2')
        #fl.show();quit()
        gridandtick(gca(),gridcolor=[1,1,1])
        fl.next('FT')
        if zero_fill:
            d.ift('t2').ft('t2',pad=2**11)
        title('FT to get $\gamma B_1/a$')
        if zero_fill:
            d.ft(ind_dim, shift=True,pad=2**11)
        else:
            d.ft(ind_dim,shift=True)
        fl.image(d[ind_dim:(-1e3*max_kHz,1e3*max_kHz)])
        #fl.show();quit()
        fl.next('absFT')
        title('FT to get $\gamma B_1/a$')
        #d.extend('t2',
        fl.image(abs(d[ind_dim:(-1e3*max_kHz,1e3*max_kHz)]))
        gridandtick(gca(), gridcolor=[1,1,1])
        fl.show();quit()
                
