from pylab import *
from pyspecdata import *
from proc_scripts import postproc_dict
zero_fill = False
with figlist_var() as fl:
    for filename,postproc,fslice,tslice,max_kHz in [
            ('201228_Ni_sol_probe_nutation_amp_2','spincore_nutation_v2',
                (-28e3,23e3),(-2e-3,2e-3),150),
            ('210120_Ni_sol_probe_nutation_amp_3','spincore_nutation_v2',
                (-16e3,15e3),(-1e-3,1e-3),200),
            ]:
        print('analyzing', filename)
        fl.basename = filename + '(preproc)\n'
        d = find_file(filename,
                exp_type='nutation',
                expno='nutation',
                postproc=postproc,
                lookup=postproc_dict,
                fl=None)
        fl.basename = filename
        plen = d.get_prop('acq_params')['p90_us']*1e-6
        if 'amp' in d.dimlabels:
            d.setaxis('amp', lambda x: x*plen)
            d.set_units('amp','s')
            ind_dim = '\\tau_p a'
            d.rename('amp',ind_dim)
        if 'p_90' in d.dimlabels:
            ind_dim = 'p_90'
        print(d.getaxis(ind_dim))
        d.extend(ind_dim,0)
        print(d.getaxis(ind_dim))
        d.ift('t2')
        print("max at",abs(d['ph1',1]['ph2',-2]).mean_all_but('t2').argmax('t2').item())
        d.setaxis('t2',lambda x: x-abs(d['ph1',1]['ph2',-2]).mean_all_but('t2').argmax('t2').item())
        fl.next('time domain -- cropped log before t-slice')
        fl.image(d.C.cropped_log(), interpolation='bilinear')
        # {{{ vlines
        oom = det_oom(d.getaxis('t2'))
        for x in tslice:
            if x is not None:
                axvline(x=x/10**oom, color='white', alpha=0.5)
        # }}}
        if tslice is not None:
            d = d['t2':tslice]
        print("signal max: %g"%abs(d['ph1',1]['ph2',-2]).data.max())
        if zero_fill:
            d.ft('t2',pad=2**11) # to improve smoothness of final result
        else:
            d.ft('t2')
        fl.next('frequency domain -- cropped log before f-slice')
        fl.image(d.C.cropped_log(), interpolation='bilinear')
        # {{{ vlines
        oom = det_oom(d.getaxis('t2'))
        for x in fslice:
            if x is not None:
                axvline(x=x/10**oom, color='white', alpha=0.5)
        # }}}
        d = d['t2':fslice]
        fl.next('frequency domain\nafter filter and t-slice')
        d = d['ph1',1]['ph2',-2]
        fl.image(d)
        gridandtick(gca(), gridcolor=[1,1,1])
        fl.next('FT')
        title('FT to get $\gamma B_1/a$')
        if zero_fill:
            d.ft(ind_dim, shift=True, pad=2**11)
        else:
            d.ft(ind_dim, shift=True)
        fl.image(d[ind_dim:(-1e3*max_kHz,1e3*max_kHz)])
        fl.next('abs FT')
        title('FT to get $\gamma B_1/a$')
        fl.image(abs(d[ind_dim:(-1e3*max_kHz,1e3*max_kHz)]))
        gridandtick(gca(), gridcolor=[1,1,1])
