from pylab import *
from pyspecdata import *
from scipy.optimize import minimize, leastsq
from sympy import exp as s_exp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, latex, Symbol
from pyspecProcScripts import *
from .simple_functions import select_pathway
from .correlation_alignment import correl_align
from .integral_w_error import integral_w_errors
from .DCCT_func import DCCT
import matplotlib.patches as patches
this_figsize = (6,12)
t2 = symbols('t2')
def process_data(s,searchstr='',
        signal_pathway={'ph1':0,'ph2':1},
        excluded_pathways = [(0,0)],
        f_range=(None,None),
        t_range=(0,0.083),
        sgn=None,
        direct='t2',
        indirect='indirect',
        Real=False,
        alias_slop=3,
        clock_correction = True,
        error_bars = True,
        correlate = True,
        fl=None):
    """Applies appropriate phase corrections, along with DC offset corrections, 
    alignment, and integrates the data to transform 2D datasets into a table
    of integrals.

    Parameters
    ==========
    searchstr:      str
                    string for the title of the plots produced.
    signal_pathway: dict
                    Dictionary of the path of the desired signal.
    excluded_pathways:  dict
                        Dictionary of all coherence pathways that are not 
                        the signal pathway.
    f_range:        tuple
                    Frequency range over which the signal resides.
    t_range:        tuple
                    Time range over which the signal resides.
    sgn:            nddata
                    A dataset with all +1 or -1 (giving the sign of the original signal)
                    Does *not* include the 'direct' dimension.
    direct:         str
                    Direct dimension of the signal
    indirect:       str
                    Indirect dimension of the signal. Usually 'power' or 'nScans' 
                    for example.
    Real:           bool
                    If true applies power axis corrections as it is real data.
    alias_slop:     int
                    Aliasing_slop used in the hermitian function.
    clock_correction:   bool
                        If true, will apply a clock correction. Is needed for IR
                        but not enhancement data.
    error_bars:     bool
                    Option to apply integration *with* error bars.
    correlate:      bool
                    Option to apply correlation alignment to data.
            
    Returns
    =======
    s:              nddata
                    Data with applied corrections but *not* integrated
    s_int:          nddata
                    Integrated and corrected data
    """                
    signal_keys = list(signal_pathway)
    signal_values = list(signal_pathway.values())
    s.ift(direct)
    if Real:
        p_axis = s.getaxis('power')
        power_axis_dBm = array(s.get_prop('meter_powers'))
        power_axis_W = zeros_like(power_axis_dBm)
        power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))
        power_axis_W = r_[0,power_axis_W]
    s.reorder([indirect,direct],first=False)
    #{{{DC offset correction
    s.ift(list(signal_pathway))
    t_start = t_range[-1] / 4
    t_start *= 3
    rx_offset_corr = s['t2':(t_start,None)]
    rx_offset_corr = rx_offset_corr.data.mean()
    s -= rx_offset_corr
    s.ft(direct)
    s.ft(list(signal_pathway),unitary=True)
    #}}}
    zero_crossing = abs(select_pathway(s[direct:f_range],signal_pathway)).C.sum(direct).argmin(indirect,raw_index=True).item()
    #{{{phase correction
    s = s[direct:f_range]
    s.ift(list(signal_pathway))
    raw_s = s.C
    s.ft(list(signal_pathway),unitary=True)
    s.ift(direct)
    s /= zeroth_order_ph(select_pathway(s,signal_pathway))
    best_shift = hermitian_function_test(select_pathway(s.C.mean(indirect)*sgn,signal_pathway),aliasing_slop=alias_slop)
    logger.info(strm("best shift is", best_shift))
    s.setaxis(direct,lambda x: x-best_shift).register_axis({direct:0})
    s.ft(direct)
    ph_corr_s = s.C
    s.ift(direct)
    if 'ph2' in s.dimlabels:
        s.reorder(['ph1','ph2',indirect,direct])
    else:
        s.reorder(['ph1',indirect,direct])
    #}}}
    #{{{Correlate
    if correlate:
        s.ft(direct)
        mysgn = determine_sign(select_pathway(s['t2':f_range], signal_pathway))
        s.ift(list(signal_pathway))
        opt_shift, sigma, my_mask = correl_align(s*mysgn,indirect_dim=indirect,
                signal_pathway=signal_pathway)
        s.ift(direct)
        s *= np.exp(-1j*2*pi*opt_shift*s.fromaxis(direct))
        s.ft(direct)
        s.ift(direct)
        s.ft(list(signal_pathway),unitary=True)
        s.ft(direct)
        if 'ph2' in s.dimlabels:
            s.reorder(['ph1','ph2',indirect,direct])
        else:
            s.reorder(['ph1',indirect,direct])
        aligned_s = s.C
        scale_factor = abs(aligned_s.C).max().item()
     #}}}  
    if fl:
        if fl.basename=='(Enhancement)':
            fl.next('')
        DCCT(raw_s,fl.next('Raw Data',figsize=(6,12)),total_spacing=0.1,
                custom_scaling=True, scaling_factor = scale_factor)
        plt.title('Raw Data for %s'%searchstr)
        DCCT(ph_corr_s,fl.next('Phased Data',figsize=(6,12)),total_spacing=0.1,
                custom_scaling=True,scaling_factor = scale_factor)
        plt.title('Phased Data for %s'%searchstr)
        if correlate:
            DCCT(aligned_s,fl.next('Aligned Data',figsize=(6,12)),total_spacing=0.1,
                    scaling_factor = scale_factor)
            plt.title('Aligned Data for %s'%searchstr)
    s.ift(direct)
    if clock_correction:
        #{{{clock correction
        clock_corr = nddata(np.linspace(-2,2,2500),'clock_corr')
        s.ft(direct)
        if 'vd' is indirect:
            s_clock = s['ph1',0]['ph2',1].sum(direct)
        else:
            s_clock = s['ph1',0].sum(direct)
        s.ift(list(signal_pathway),unitary=True)
        min_index = abs(s_clock).argmin(indirect,raw_index=True).item()
        s_clock *= np.exp(-1j*clock_corr*s.fromaxis(indirect))
        s_clock[indirect,:min_index+1] *= -1
        s_clock.sum(indirect).run(abs)
        clock_corr = s_clock.argmax('clock_corr').item()
        s *= np.exp(-1j*clock_corr*s.fromaxis(indirect))
        s.ft(list(signal_pathway),unitary=True)
        if fl:
            DCCT(s,fl.next('After Auto-Clock Correction'),total_spacing=0.2,
                    custom_scaling=True, scaling_factor = scale_factor)
            plt.title('After Auto-Clock Correction %s'%searchstr)
        s.ift(direct)   
        #}}}
    s_after = s.C 
    s_after[indirect,zero_crossing] *= -1
    s_after.ft(direct)
    s_after.ift(direct)
    s_after = s_after[direct:(0,None)]
    s_after[direct,0] *= 0.5
    s_after.ft(direct)
    if fl:
        DCCT(s_after,fl.next('FID',figsize=this_figsize), total_spacing=0.2,
                scaling_factor = scale_factor)
        plt.title('FID')
    if 'ph2' in s.dimlabels:
        error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
                - set(excluded_pathways)
                - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
        error_path = [{'ph1':j,'ph2':k} for j,k in error_path]
    else:
        error_path = (set(((j) for j in range(ndshape(s)['ph1'])))
                - set(excluded_pathways)
                - set([(signal_pathway['ph1'])]))
        error_path = [{'ph1':j} for j in error_path]
            #{{{Integrate with error bars
    if error_bars:
        s_int,frq_slice = integral_w_errors(s_after,signal_pathway,error_path,
                convolve_method='Gaussian',
                indirect=indirect, return_frq_slice = True)
        x = s_int.get_error()
        x[:] /= sqrt(2)
        if fl is not None:
            fl.next('Real with Integration Bounds',figsize=this_figsize)
            left_pad, bottom_pad, width_pad, top_pad = DCCT(s,fl.next('Real with Integration Bounds',figsize=this_figsize),just_2D=True)
            x = s_after.getaxis(direct)
            dx = x[1]-x[0]
            y = s_after.getaxis(indirect)
            dy = y[1]-y[0]
            start_y = s_after.getaxis(indirect)[0]
            stop_y = s_after.getaxis(indirect)[-1]
            #Does not plot properly if below figure is not produced
            ax = plt.axes([left_pad,bottom_pad,width_pad,top_pad])
            yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,10])
            majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=2, steps=[1,10])
            minorLocator = lambda: mticker.AutoMinorLocator(n=4)
            ax.xaxis.set_major_locator(majorLocator())
            ax.xaxis.set_minor_locator(minorLocator())
            ax.set_ylabel(None)
            fl.image(select_pathway(s_after.real.run(complex128),signal_pathway),
                black=False)
            x1 = x[0]
            y1 = start_y - dy
            x2 = frq_slice[-1]
            tall = (y[-1]+dy)-(y[0]-dy)
            wide = frq_slice[0] - x1 - dx
            wide2 = (x[-1]+dx)-x2
            p = patches.Rectangle((x1,y1),wide,tall,angle=0.0,linewidth=1,fill=None,
                    hatch='//',ec='k')
            ax.add_patch(p)
            q = patches.Rectangle((x2,y1),wide2,tall,angle=0.0,linewidth=1,fill=None,
                    hatch='//',ec='k')
            ax.add_patch(q)
    if indirect is 'vd':
        s_int = s_int
        if fl is not None:
            fig2=figure(figsize=this_figsize)
            fl.next('Integrated Data', fig=fig2)
            fl.plot(s_int,'o')
    else:
        s_int[indirect,:] /= s_int.data[0]
        if Real:
            s_int.setaxis('power',power_axis_W)
        if fl is not None:
            fig1 = figure(figsize=this_figsize)
            fl.next('Integrated Data',fig=fig1)
            fl.plot(s_int['power',:-3],'o')
    #}}}    
    return s_int, s
    
