from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from sympy import symbols
import matplotlib.pyplot as plt
fl = fl_mod()
t2 = symbols('t2')
filter_bandwidth = 20e3
gamma_eff = (14.904100/3507.57) # MHz / G
for thisfile,exp_type,nodename,postproc,label_str,freq_slice in [
        ('220126_150mM_TEMPOL_field_dep','ODNP_NMR_comp/field_dependent',
            'field_sweep','field_sweep_v1',
            'TEMPO field sweep',(-250,250)),
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup = lookup_table)
    #{{{Obtain ESR frequency and chunk/reorder dimensions
    nu_B12 = (s.get_prop('acq_params')['mw_freqs'][0])/1e9
    #}}}
    #{{{DC offset correction
    t2_max = s.getaxis('t2')[-1]
    rx_offset_corr = s['t2':(t2_max*0.75,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft(['ph1'])
    fl.next('raw data -- coherence channels')
    print(ndshape(s))
    s.reorder(['ph1','indirect','power','t2'])
    fl.image(s.C.setaxis(
'indirect','#').set_units('indirect','scan #'))
    #}}}
    #{{{frequency filtering and rough center
    s.ft('t2',shift=True)
    fl.next('frequency domain raw data')
    fl.image(s.C.setaxis(
'indirect','#').set_units('indirect','scan #'))
    s = s['t2':(-1e3,1e3)]
    s=s['ph1',1]['power',0]['nScans',0].C
    s.ift('t2')
    s.ft('t2')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    best_shift = hermitian_function_test(s)
    s.setaxis('t2',lambda x: x-best_shift).register_axis({'t2':0})
    s.ft('t2')
    fl.next('line plots')
    nu_NMR = []
    offsets = []
    s = s['indirect',:-1]
    for z in range(len(s.getaxis('indirect')[:]['Field'])):
        fl.plot(abs(s['indirect',z]),label = '%0.2f'%s.getaxis('indirect')[z]['Field'])
        offset = abs(s['indirect',z].C).argmax('t2')
        offsets.append(offset)
        true_carrier_freq = s.getaxis('indirect')[z]['carrierFreq']
        nu_rf = true_carrier_freq*1e6 + (offset)
        nu_rf /= 1e6
        plt.axvline(x=freq_slice[0])
        plt.axvline(x=freq_slice[-1])
        nu_NMR.append(nu_rf.data) 
    s = s['t2':freq_slice].sum('t2')
    #}}}
    #{{{convert x axis to ppt = nu_NMR/nu_ESR
    ppt = (nu_NMR / nu_B12)
    s.setaxis('indirect',ppt)
    s.rename('indirect','ppt')
    #}}}
    #{{{Fitting
    fl.next('sweep, without hermitian')
    fl.plot(abs(s),'-o')
    field_idx = (abs(s.data)).argmax()
    fitting = abs(s).polyfit('ppt',order=2)
    x_min = s.getaxis('ppt')[0]
    x_max = s.getaxis('ppt')[-1]
    Field = nddata(r_[x_min:x_max:100j],'ppt')
    fl.plot(Field.eval_poly(fitting,'ppt'),label='fit')
    #}}}
    logger.info(strm("ESR frequency is %f"%(nu_B12)))
    logger.info(strm('The fit finds a max with ppt value:',
        Field.eval_poly(fitting,'ppt').argmax().item()))
    logger.info(strm('The data finds a ppt value', abs(s).argmax().item()))
fl.show()
