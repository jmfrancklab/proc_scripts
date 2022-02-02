"""Check ppt value using a field sweep
======================================

A DC offset correction, initial frequency slice 
and rough centering are applied to a previously
acquired TEMPOL dataset. It is then sliced more based
on the 1D diagnostic plot labeled 'Field Slicing'. The 
field is converted to ppt using the frequency stored as
an acquisition parameter in the dataset and the effective
gamma used that day. This will verify the previously required
ppt value for a given spin label.
"""

from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from sympy import symbols
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 20e3
gamma_eff = (14.904100/3507.57) # MHz / G
rcParams["image.aspect"] = "auto" # needed for sphinx gallery

#sphinx_gallery_thumbnail_number = 1


for thisfile,exp_type,nodename,postproc,label_str,freq_slice in [
        ('220126_150mM_TEMPOL_field_dep','ODNP_NMR_comp/field_dependent',
            'field_sweep','field_sweep_v1',
            'TEMPOL field sweep',(-250,250)),
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=lookup_table)
    fig, ax_list = subplots(1, 3)
    #{{{Obtain ESR frequency and chunk/reorder dimensions
    nu_B12 = (s.get_prop('acq_params')['mw_freqs'][0])/1e9
    #}}}
    #{{{DC offset correction
    t2_max = s.getaxis('t2')[-1]
    rx_offset_corr = s['t2':(t2_max*0.75,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft(['ph1'])
    s.reorder(['ph1','indirect','power','t2'])
    fl.next('Field Sweep Processing',fig=fig)
    #}}}
    #{{{frequency filtering and rough center
    s.ft('t2',shift=True)
    fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'),ax = ax_list[0])
    ax_list[0].set_title('Raw data\nFrequency Domain')
    s = s['t2':(-1e3,1e3)]
    s=s['ph1',1]['power',0]['nScans',0].C
    s.ift('t2')
    s.ft('t2')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    best_shift = hermitian_function_test(s)
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    s.ft('t2')
    nu_NMR=[]
    offsets = []
    s = s['indirect',:-1]
    for z in range(len(s.getaxis('indirect')[:]['Field'])):
        fl.plot(abs(s['indirect',z]),label = '%0.2f'%s.getaxis('indirect')[z]['Field'],
                ax=ax_list[1])
        offset = abs(s['indirect',z].C).argmax('t2')
        offsets.append(offset)
        true_carrier_freq = s.getaxis('indirect')[z]['carrierFreq']
        nu_rf = true_carrier_freq*1e6 + (offset)
        nu_rf /= 1e6
        nu_NMR.append(nu_rf.data) 
        ax_list[1].axvline(x = freq_slice[0])
        ax_list[1].axvline(x=freq_slice[-1])
    ax_list[1].set_title('Field Slicing')
    s = s['t2':freq_slice].sum('t2')
    #}}}
    #{{{convert x axis to ppt = v_NMR/v_ESR
    ppt = nu_NMR / nu_B12 
    s.setaxis('indirect',ppt)
    s.rename('indirect','ppt')
    #}}}
    #{{{Fitting
    fl.plot(abs(s),'o-',ax=ax_list[2])
    ax_list[2].set_title('Sweep, \nwithout a Hermitian')
    field_idx = (abs(s.data)).argmax()
    fitting = abs(s).polyfit('ppt',order=2)
    x_min = s.getaxis('ppt')[0]
    x_max = s.getaxis('ppt')[-1]
    Field = nddata(r_[x_min:x_max:100j],'ppt')
    fl.plot(Field.eval_poly(fitting,'ppt'),label='fit',ax=ax_list[2])
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    #}}}
    logger.info(strm("ESR frequency is %f"%(nu_B12)))
    logger.info(strm('The fit finds a max with ppt value:',
        Field.eval_poly(fitting,'ppt').argmax().item()))
    logger.info(strm('The data finds a ppt value', abs(s).argmax().item()))
fl.show()
