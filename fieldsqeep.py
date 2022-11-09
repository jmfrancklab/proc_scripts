from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from sympy import symbols
fl = figlist_var()
signal_pathway = {'ph1':1}
for thisfile,exp_type,nodename,postproc,label_str in [
        ('220818_pRA174_MSL_500uM_field.h5',
            'ODNP_NMR_comp/field_dependent',
            'field_1',
            'field_sweep_v2',
            'TEMPOL field sweep')
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc, lookup=lookup_table)
    #{{{Obtain ESR frequency and chunk/reorder dimensions
    nu_B12 = s.get_prop('acq_params')['mw_freqs']/1e9
    s.reorder(['ph1','indirect','nScans','t2'])
    #}}}
    #{{{DC offset correction
    s.ift('t2')
    s.ift(['ph1'])
    t2_max = s.getaxis('t2')[-1]
    rx_offset_corr = s['t2':(t2_max*0.75,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    fl.next('time domain')
    s.ft(['ph1'])
    fl.image(s['t2':(None,0.05)])
    s.ft('t2')
    #}}}
    fl.next('frequency domain')
    fl.image(s)
    #{{{frequency filtering and rough center
    s = s['t2':(-1.0e3,1.5e3)]
    fl.next('before phasing')
    fl.plot(s.C.mean('nScans'))
    s.ift('t2')
    s.set_units('t2','s')
    best_shift = 0.0035
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    s /= zeroth_order_ph(select_pathway(s['t2',0],signal_pathway))
    s = select_pathway(s,signal_pathway)
    s.ft('t2')
    fl.next('phased')
    fl.plot(s.C.mean('nScans'))
    fl.next('phased im')
    fl.image(s)
    nu_NMR=[]
    all_offsets = zeros(len(s.getaxis('indirect')))
    for z in range(len(s.getaxis('indirect'))):
        fl.next('Field slicing')
        fl.plot(s['indirect',z],label='scan %d'%z)
        offset = s['indirect',z].C.mean('nScans').argmax('t2').item()
        all_offsets[z] = offset
        carrier_freq_MHz = s.getaxis('indirect')[z]['carrierFreq']
        print("for indirect %d"%z)
        print(" your offset is %.4f"%offset)
        nu_rf = carrier_freq_MHz - offset/1e6
        nu_NMR.append(nu_rf) 
    s.ift('t2')
    s *= exp(-1j*2*pi*nddata(all_offsets, [-1],['indirect'])*s.fromaxis('t2'))
    s.ft('t2')
    frq_slice = s.C.mean('nScans').mean('indirect').contiguous(lambda x: x.real > 0.05*s.real.data.max())[0]
    fl.next('phased data')
    fl.plot(s.C.mean('nScans'))
    plt.axvline(x = frq_slice[0])
    plt.axvline(x = frq_slice[-1])
    s = s['t2':frq_slice].mean('nScans').integrate('t2')
    #}}}
    #{{{convert x axis to ppt = v_NMR/v_ESR
    ppt = np.asarray(nu_NMR) / np.asarray(nu_B12)
    ppt.sort()
    s.setaxis('indirect',ppt)
    #s.set_units('indirect','MHz/GHz')
    s.rename('indirect','ppt')
    #}}}
    #{{{Fitting
    fl.next('fit')
    fl.plot(s,'o')
    fitting = s.polyfit('ppt',order=4)
    x_min = s['ppt'][0]
    x_max = s['ppt'][-1]
    Field = nddata(r_[x_min:x_max:100j],'ppt')
    fl.plot(Field.eval_poly(fitting,'ppt'),label='fit')
    #}}}
    print("ESR frequency is %f"%(nu_B12))
    print('The fit finds a max with ppt value:',
        Field.eval_poly(fitting,'ppt').argmax().item())
    print('The data finds a ppt value', abs(s).argmax().item())
fl.show()
