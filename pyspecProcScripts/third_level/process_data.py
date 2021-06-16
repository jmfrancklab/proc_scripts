from pyspecdata import *
from pyspecProcScripts import *
def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
def proc_data(s,label='',indirect = 'vd',fl=None,signal_pathway={'ph1':0,'ph2':1},
        excluded_pathways = [(0,0)], clock_correction=False,flip=False,
        f_range=(None,None), t_range=(None,83e-3),sign=None):
    s *= sign
    if 'ph2' in s.dimlabels:
        s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
        s.ift('t2')
        #{{{ Applying DC offset
        s.ift(['ph1','ph2'])
        t_start = t_range[-1] / 4
        t_start *= 3
        rx_offset_corr = s['t2':(t_start,None)]
        rx_offset_corr = rx_offset_corr.data.mean()
        s -= rx_offset_corr
        s.ft('t2')
        s.ft(['ph1','ph2'])
        #}}}
        s = s['t2':f_range]
        s.ift('t2')
        if clock_correction:
            #{{{ clock correction
            clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
            s.ft('t2')
            if fl is not None:
                fl.next('before clock correction')
                fl.image(s.C.setaxis(indirect,'#').set_units(indirect,'scan #'))
            s_clock=s['ph1',1]['ph2',0].sum('t2')
            s.ift(['ph1','ph2'])
            min_index = abs(s_clock).argmin(indirect,raw_index=True).item()
            s_clock *= np.exp(-1j*clock_corr*s.fromaxis('vd'))
            s_clock[indirect,:min_index+1] *=-1
            s_clock.sum(indirect).run(abs)
            if fl is not None:
                fl.next('clock correction')
                fl.plot(s_clock,'.',alpha=0.7)
            clock_corr = s_clock.argmax('clock_corr').item()
            plt.axvline(x=clock_corr, alpha=0.5, color='r')
            s *= np.exp(-1j*clock_corr*s.fromaxis(indirect))
            s.ft(['ph1','ph2'])
            if fl is not None:
                fl.next('after auto-clock correction')
                fl.image(s.C.setaxis(indirect,'#'))
            s.ift('t2')
        #{{{Applying phase corrections    
        best_shift,max_shift = hermitian_function_test(select_pathway(s.C.mean(indirect),signal_pathway))
        logger.info(strm("best shift is", best_shift))
        s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
        if fl is not None:
            fl.next('time domain after hermitian test')
            fl.image(s.C.setaxis(indirect,'#').set_units(indirect,'scan #'))
        s.ft('t2')
        if fl is not None:
            fl.next('frequency domain after hermitian test')
            fl.image(s.C.setaxis(indirect,'#').set_units(indirect,'scan #'))
            #}}}
        s.ift('t2')
        s.ift(['ph1','ph2'])
        phasing = s['t2',0].C
        phasing.data *= 0
        phasing.ft(['ph1','ph2'])
        phasing['ph1',0]['ph2',1] = 1
        phasing.ift(['ph1','ph2'])
        ph0 = s['t2':0]
        ph0 /= abs(ph0)
        s /= ph0
        s.ft(['ph1','ph2'])
        s.ft('t2')
        s.reorder(['ph1','ph2',indirect,'t2'])
        #{{{Correlation Alignment
        fl.basename='correlation subroutine:'
        #for the following, should be modified so we can pass a mask, rather than specifying ph1 and ph2, as here
        s_aligned,opt_shift,sigma = correl_align(s,ODNP=False,indirect_dim=indirect,
                ph1_selection=signal_pathway['ph1'],ph2_selection=signal_pathway['ph2'],
                sigma=10)
        s = s_aligned
        fl.basename = None
        if fl is not None:
            fl.next(r'after correlation, $\varphi$ domain')
            fl.image(s.C.setaxis(indirect,'#').set_units(indirect,'scan #'))   
        s.ift('t2')
        s.ft(['ph1','ph2'])
        if fl is not None:
            fl.next(r'after correlation')
            fl.image(s.C.setaxis(indirect,'#').set_units(indirect,'scan #'))
        if 'nScans' in s.dimlabels:
            s.mean('nScans')  
        s.ft('t2')
        s.ift('t2')
        s = s['t2':(0,t_range[-1])]
        s['t2':0] *= 0.5
        s.ft('t2')
        zero_crossing=abs(select_pathway(s,signal_pathway)).sum('t2').argmin(indirect,raw_index=True).item()
        if flip:
            s[indirect,:zero_crossing] *= -1
        # {{{ this is the general way to do it for 2 pulses I don't offhand know a compact method for N pulses
        error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
                - set(excluded_pathways)
                - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
        error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
        # }}}
        #{{{Integrating with associated error from excluded pathways    
        s_int,frq_slice = integral_w_errors(s,signal_pathway,error_path,
                fl=fl,return_frq_slice=True)
        x = s_int.get_error()
        x[:] /= sqrt(2)
        return s_int

    else:
        s *= sign
        s['ph1',0]['t2':0] = 0
        #{{{ Applying DC offset
        s.ift(['ph1'])
        t_start = t_range[-1] / 4
        t_start *= 3
        rx_offset_corr = s['t2':(t_start,None)]
        rx_offset_corr = rx_offset_corr.data.mean()
        s -= rx_offset_corr
        s.ft('t2')
        s.ft(['ph1'])
        if 'nScans' in s.dimlabels:
            s.reorder(['ph1','nScans',indirect,'t2'])
        else:
            s.reorder(['ph1',indirect,'t2'])
        #}}}
        s = s['t2':f_range]
        s.ift('t2')
        if clock_correction:
            #{{{ clock correction
            clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
            s.ft('t2')
            if fl is not None:
                fl.next('before clock correction')
                fl.image(s)
            s_clock=s['ph1',1].sum('t2')
            s.ift(['ph1'])
            min_index = abs(s_clock).argmin(indirect,raw_index=True).item()
            s_clock *= np.exp(-1j*clock_corr*s.fromaxis(indirect))
            s_clock[indirect,:min_index+1] *=-1
            s_clock.sum(indirect).run(abs)
            if fl is not None:
                fl.next('clock correction')
                fl.plot(s_clock,'.',alpha=0.7)
            clock_corr = s_clock.argmax('clock_corr').item()
            plt.axvline(x=clock_corr, alpha=0.5, color='r')
            s *= np.exp(-1j*clock_corr*s.fromaxis(indirect))
            s.ft(['ph1'])
            if fl is not None:
                fl.next('after auto-clock correction')
                fl.image(s.C.setaxis(indirect,'#'))
            s.ift('t2')
        #{{{Applying phase corrections    
        best_shift,max_shift = hermitian_function_test(select_pathway(s,signal_pathway).C.convolve('t2',0.001))
        logger.info(strm("best shift is", best_shift))
        s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
        if fl is not None:
            fl.next('time domain after hermitian test')
            fl.image(s.C.setaxis(
'nScans','#').set_units('nScans','scan #').setaxis(
'power','#').set_units('power','scan #'))
        s.ft('t2')
        if fl is not None:
            fl.next('frequency domain after hermitian test')
            fl.image(s.C.setaxis(
'nScans','#').set_units('nScans','scan #').setaxis(
'power','#').set_units('power','scan #'))
            #}}}
        s.ift('t2')
        s.ift(['ph1'])
        phasing = s['t2',0].C
        phasing.data *= 0
        phasing.ft(['ph1'])
        phasing['ph1',1] = 1
        phasing.ift(['ph1'])
        ph0 = s['t2':0]/phasing
        ph0 /= abs(ph0)
        s /= ph0
        s.ft(['ph1'])
        s.ft('t2')
        if 'nScans' in s.dimlabels:
            s.reorder(['ph1','nScans',indirect,'t2'])
        else:
            s.reorder(['ph1',indirect,'t2'])
        #{{{Correlation Alignment
        fl.basename='correlation subroutine:'
        #for the following, should be modified so we can pass a mask, rather than specifying ph1 and ph2, as here
        if 'nScans' in s.dimlabels:
            logging.info(strm("Aligning per transient"))
            s_aligned,opt_shift,sigma = correl_align(s,indirect_dim='nScans',
                    ph1_selection=signal_pathway['ph1'],
                    sigma=0.001)
        else:
            s_aligned,opt_shift,sigma = correl_align(s,indirect_dim=indirect,
                ph1_selection=1,sigma=0.001)
        s = s_aligned
        fl.basename = None
        if fl is not None:
            fl.next(r'after correlation, $\varphi$ domain')
            fl.image(s.C.setaxis(
'nScans','#').set_units('nScans','scan #').setaxis(
'power','#').set_units('power','scan #'))   
        s.ift('t2')
        s.ft(['ph1'])
        if 'nScans' in s.dimlabels:
            s.mean('nScans')
        if fl is not None:
            fl.next(r'after correlation')
            fl.image(s.C.setaxis(
'power','#').set_units('power','scan #'))  
        s.ft('t2')
        s.ift('t2')
        s = s['t2':(0,t_range[-1])]
        s['t2':0] *= 0.5
        s.ft('t2')
        error_pathway = (set(((j) for j in range(ndshape(s)['ph1'])))
                - set(excluded_pathways)
                - set([(signal_pathway['ph1'])]))
        error_pathway = [{'ph1':j} for j in error_pathway]
        s_,frq_slice = integral_w_errors(s,signal_pathway,error_pathway,
                indirect='power', fl=fl, return_frq_slice=True)
        x = s_.get_error()
        x[:] /= sqrt(2)
        s = s_.C
    return s
  
   

    

       
        



