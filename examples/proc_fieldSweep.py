"""Check NMR/ESR resonance ratio using a field sweep
====================================================

Analyzes field sweep data.  Determines the resonance frequency from the carrier
frequency and the offset of the signal, converts to MHz, then divide by the
Bridge 12 μw frequency stored in the file to get the resonance ratio of MHz/GHz
"""

from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from sympy import symbols
fl = figlist_var()
rcParams["image.aspect"] = "auto" # needed for sphinx gallery

#sphinx_gallery_thumbnail_number = 1

signal_pathway = {'ph1':1}
for thisfile,exp_type,nodename,postproc,label_str,freq_slice in [
        ('220217_5mM_TEMPOL_field_dep',
            'ODNP_NMR_comp/field_dependent',
            'Field_sweep',
            'field_sweep_v1',# in newer version files, we should be setting
            #                  postproc_type as a top-level attribute of the
            #                  nddata, and this should not be necessary
            'TEMPOL field sweep',
            (-250,250)),
        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=lookup_table)
    #{{{Obtain ESR frequency and chunk/reorder dimensions
    nu_B12_GHz = (s.get_prop('acq_params')['mw_freqs'][0])/1e9
    #}}}
    #{{{DC offset correction
    s.ift('t2')
    s.ift('ph1')
    t2_max = s.getaxis('t2')[-1]
    rx_offset_corr = s['t2':(t2_max*0.75,None)]
    rx_offset_corr = rx_offset_corr.mean(['t2'])
    s -= rx_offset_corr
    s.ft(['ph1'])
    s.ft('t2')
    #}}}
    # {{{ set up figure and plot raw data
    fig, ax_list = subplots(1, 3)
    fl.next('Field Sweep Processing',fig=fig)
    fl.image(s,ax = ax_list[0])
    ax_list[0].set_title('Raw data\nFrequency Domain')
    # }}}
    #{{{frequency filtering and phase correct
    s = s['t2':(-1e3,1e3)]
    s.ift('t2')
    best_shift = hermitian_function_test(select_pathway(s,signal_pathway))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    s.ft('t2')
    s = select_pathway(s,signal_pathway)
    s /= zeroth_order_ph(s.C.mean('t2'))
    ## JF review to here
    nu_NMR=[]
    for z in range(len(s.getaxis('indirect')[:]['Field'])):
        fl.plot(abs(s['indirect',z]),ax=ax_list[1]) #there is some rolling of baseline without abs
        offset = abs(s['indirect',z].C.mean('nScans')).argmax('t2')
        true_carrier_freq = s.getaxis('indirect')[z]['carrierFreq']
        nu_rf = true_carrier_freq*1e6 + (offset) #in Hz
        nu_rf /= 1e6 #back to MHz
        nu_NMR.append(nu_rf.data) 
        ax_list[1].axvline(x = freq_slice[0])
        ax_list[1].axvline(x=freq_slice[-1])
    ax_list[1].set_title('Field Slicing')
    s = s['t2':freq_slice].mean('nScans').sum('t2')
    #}}}
    #{{{convert x axis to ppt = v_NMR/v_ESR
    ppt = nu_NMR / nu_B12_GHz 
    s.setaxis('indirect',ppt)
    s.rename('indirect','ppt')
    #}}}
    #{{{Fitting
    fl.plot(s,'o-',ax=ax_list[2])
    ax_list[2].set_title('Field Sweep ppt')
    fitting = s.polyfit('ppt',order=2)
    Field = nddata(r_[s.getaxis('ppt')[0]:s.getaxis('ppt')[-1]:100j],'ppt')
    fl.plot(Field.eval_poly(fitting,'ppt'),label='fit',ax=ax_list[2])
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    #}}}
    logger.info(strm("ESR frequency is %f"%(nu_B12_GHz)))
    logger.info(strm('The fit finds a max with ppt value:',
        Field.eval_poly(fitting,'ppt').argmax().item()))
    logger.info(strm('The data finds a ppt value', s.argmax().item()))
fl.show()
