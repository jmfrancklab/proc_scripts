from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from sympy import symbols
import matplotlib.pyplot as plt
fl = fl_mod()
t2 = symbols('t2')
filter_bandwidth = 20e3
gamma_eff = (14.903537/3507.48) # MHz / G
for thisfile,exp_type,nodename,postproc,label_str,freq_slice in [
        ('220112_150uM_TEMPOL_field_dep_1','ODNP_NMR_comp/field_dependent',
            'field_sweep_1','field_sweep_v1',
            'TEMPO field sweep',(-250,250)),
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup = lookup_table)
    #{{{Obtain ESR frequency and chunk/reorder dimensions
    v_B12 = (s.get_prop('acq_params')['mw_freqs'][0])/1e9
    #}}}
    #{{{DC offset correction
    t2_max = s.getaxis('t2')[-1]
    rx_offset_corr = s['t2':(t2_max*0.75,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft(['ph1'])
    fl.next('raw data -- coherence channels')
    s.reorder(['ph1','Field','power','t2'])
    fl.image(s)
    #}}}
    #{{{frequency filtering and rough center
    s.ft('t2',shift=True)
    fl.next('frequency domain raw data')
    fl.image(s)
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
    v_NMR = []
    offsets = []
    for z in range(len(s.getaxis('Field'))):
        fl.plot(abs(s['Field',z]),label = '%0.2f'%s.getaxis('Field')[z])
        offset = abs(s['Field',z].C).argmax('t2')
        offsets.append(offset)
        true_carrier_freq = s.getaxis('Field')[z] * gamma_eff
        v_rf = (s.get_prop('acq_params')['carrierFreq_MHz']*1e6 + (offset))/1e6
        #plt.axvline(x=freq_slice[0])
        #plt.axvline(x=freq_slice[-1])
        v_NMR.append(v_rf.data) 
    s_ = s['t2':freq_slice].sum('t2')
    print("s.getaxis('Field') =", s.getaxis('Field'))
    print("offsets:",offsets)
    #fl.show();quit()
    #}}}
    #{{{convert x axis to ppt = v_NMR/v_ESR
    ppt = (v_NMR / v_B12)
    s_.setaxis('Field',ppt)
    s_.rename('Field','ppt')
    #}}}
    #{{{Fitting
    fl.next('sweep, without hermitian')
    fl.plot(abs(s_),'-o')
    fl.show();quit()
    field_idx = (abs(s_.data)).argmax()
    fitting = abs(s_).polyfit('ppt',order=2)
    x_min = s_.getaxis('ppt')[0]
    x_max = s_.getaxis('ppt')[-1]
    Field = nddata(r_[x_min:x_max:100j],'ppt')
    fl.plot(Field.eval_poly(fitting,'ppt'),label='fit')
    #}}}
    logger.info(strm("ESR frequency is %f"%(v_B12)))
    logger.info(strm('The fit finds a max with ppt value:',
        Field.eval_poly(fitting,'ppt').argmax().item()))
    logger.info(strm('The data finds a ppt value', abs(s_).argmax().item()))
fl.show()
