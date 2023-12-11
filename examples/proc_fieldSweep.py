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
init_logging(level='debug')
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
    fig, ax_list = subplots(1, 4, figsize=(10,3))
    fl.next('Field Sweep Processing',fig=fig)
    fl.image(s,ax = ax_list[0])
    ax_list[0].set_title('Raw data\nFrequency Domain')
    # }}}
    #{{{frequency filtering and phase correct
    s = s['t2':(-1e3,1e3)]
    s.ift('t2')
    s.set_units('t2','s')
    best_shift,_ = hermitian_function_test(select_pathway(s,signal_pathway),
            echo_before = s.get_prop('acq_params')['tau_us']*1e-6*1.5)
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    s = select_pathway(s,signal_pathway)
    s /= zeroth_order_ph(s.C.mean('t2'))
    nu_NMR=[]
    assert set(s.getaxis('indirect').dtype.names) == {'Field', 'carrierFreq'}, "'indirect' axis should be a structured array that stores the carrier frequency and the field"
    # {{{ the following works fine for all the cases that we've tested, and so
    # there's no need to change it.  But, after resolving everything, realized
    # that the best way to deal with this would be to run a correlation on all
    # the scans to determine the relative offsets, which would give both the
    # relative offsets and the frequency shifts needed to align the signals, at
    # the end.
    all_offsets = zeros(len(s.getaxis('indirect')))
    for z in range(len(s.getaxis('indirect'))):
        fl.plot(s['indirect',z],ax=ax_list[1],human_units = False)
        if z == 0:
            ax_list[1].axvline(color='k',x = freq_slice[0])
            ax_list[1].axvline(color='k',x=freq_slice[-1])
        offset = s['indirect',z].C.mean('nScans').argmax('t2').item()
        all_offsets[z] = offset
        carrier_freq_MHz = s.getaxis('indirect')[z]['carrierFreq']
        field = s.getaxis('indirect')[z]['Field']
        nu_rf = carrier_freq_MHz - offset/1e6 # because in order to make my
        #                                       offset more negative, I
        #                                       increase my field, or decrease
        #                                       my carrier -- in other words
        #                                       the frequency axis for our FID
        #                                       is the negative of (ν_RF-ν_carrier)
        nu_NMR.append(nu_rf) 
    ax_list[1].set_title('Field Slicing')
    s.ift('t2')
    s *= exp(-1j*2*pi*nddata(all_offsets, [-1], ['indirect'])*s.fromaxis('t2'))
    s.ft('t2')
    # }}}
    fl.plot(s.C.mean('nScans'), ax = ax_list[2],human_units = False)
    s = s['t2':freq_slice].mean('nScans').integrate('t2')
    #}}}
    #{{{convert x axis to ppt = v_NMR/v_ESR
    ppt = nu_NMR / nu_B12_GHz 
    s.setaxis('indirect',ppt)
    s.set_units('indirect','MHz/GHz')
    s.rename('indirect','resonance ratio')
    #}}}
    #{{{Fitting
    fl.plot(s,'o-',ax=ax_list[3])
    ax_list[3].set_title('Field Sweep ppt')
    fitting = s.polyfit('resonance ratio',order=2)
    field_fine = nddata(r_[s.getaxis('resonance ratio')[0]:s.getaxis('resonance ratio')[-1]:100j],'resonance ratio')
    polyline = fl.plot(field_fine.eval_poly(fitting,'resonance ratio'),label='fit',ax=ax_list[3])
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    _, B, A = fitting
    center = -B/2/A
    #}}}
    logger.info(strm("ESR frequency is %f"%(nu_B12_GHz)))
    logger.info(strm("MHz/GHz value at peak/dip appears to be:",center))
    ax_list[3].axvline(x=center, color=polyline[0].get_color(), ls=':')
fl.show()
