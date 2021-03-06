from pyspecdata import *
from scipy.optimize import minimize,leastsq
from proc_scripts import *
from proc_scripts import postproc_dict
do_slice = False # slice frequencies and downsample -- in my hands, seems to decrease the quality of the fit 
standard_cost = False # use the abs real to determine t=0 for the blip -- this actually doesn't seem to work, so just use the max
show_transfer_func = False # show the transfer function -- will be especially useful for processing non-square shapes
logger = init_logging('info')
#init_logging(level='debug')
# 2to3 JF 1/31

fl = figlist_var()
 # {{{ load data, set units, show raw data
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        #('180806','pulse_reflection',True),
        #('181001','sprobe_t2',True),
        #('181001','sprobe_t4',True),
        #('181103','probe',True),
        #('200110','pulse_2',True),
        #('200312','chirp_coile_4',True),
        ('200103_pulse_1','test_equip','capture1','square_wave_capture_v1',True),
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 

    fl.next('Raw signal %s'%searchstr)
    fl.plot(d['ch',0], alpha=0.5, label='control') # turning off human units forces plot in just V
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    # }}}
    
    # {{{ determining center frequency and convert to
    # analytic signal, show analytic signal
    d.ft('t',shift=True) #Fourier Transform into freq domain
    d = 2*d['t':(0,None)] # throw out negative frequencies and low-pass (analytic signal -- 2 to preserve magnitude)
    #to negated the "1/2" in "(1/2)(aexp[iwt]+a*exp[-iwt])
    # ALSO -- everything should be less than 40 MHz for sure
    if do_slice:
        d = d['t':(0,40e6)]
    else:
        d['t':(0,3e6)] = 0 # filter out low-frq noise, which becomes high-frq noise on demod
    center_frq = abs(d['ch',0]).argmax('t').item() # the center frequency is now the max of the freq peak    
    logger.info(strm(("initial guess at center frequency at %0.5f MHz"%(center_frq/1e6))))
    logger.info(strm(center_frq))
    fl.next('frequency domain\n%s'%searchstr)
    fl.plot(abs(d['t':(None,40e6)]),label='Raw signal in freq domain,\nshows a bandwidth of about 20 MHz', alpha=0.5)
    axvline(x=center_frq/1e6)
    d.setaxis('t',lambda x: x-center_frq).register_axis({'t':0})
    if do_slice:
        d = d['t':(-10e6,10e6)]
    d.ift('t') #Inverse Fourier Transform back to time domain to display the decaying exponential
    fl.next('Absolute value of analytic signal, %s'%searchstr)
    fl.plot(abs(d['ch',0]), alpha=0.5, label='control') #plot the 'envelope' of the control 
    fl.plot(abs(d['ch',1]), alpha=0.5, label='reflection') #plot the 'envelope' of the reflection so no more oscillating signal
    # }}}
    
    #{{{determine the start and stop points for both the pulse, as well as the two tuning blips
    pulse_range = abs(d['ch',0]).contiguous(lambda x:  # returns list of limits for which the lambda function holds true
            x > 0.5*x.data.max())                      # So will define pulse_range as all x values where the signal goes 
                                                       # above half the max of the signal
    #}}}

    #{{{ filter for ranges >0.1 μs -- use the compact list comprehension
    def filter_range(x): return array([j for j in x if
        diff(j).item() > 0.1e-6])
    pulse_range = filter_range(pulse_range)
    if not pulse_range.shape[0] == 1:
        logger.info(strm(("seems to be more than one pulse -- on starting at " 
                + ','.join(('start '+str(j[0])+' length '+str(diff(j)) for j in pulse_range)))))   # If there is more than one section that goes above half max
    # it assumes theres more than one pulse 
    pulse_range = pulse_range[0,:]
    #}}}

    #{{{plotting reflection blip
    fl.plot(abs(d['ch',0]['t':tuple(pulse_range)]), alpha=0.1, color='k',  #shades in the section of pulse range (above half max) for 
            linewidth=10)                                                  #control 
    refl_blip_ranges = abs(d['ch',1]).contiguous(lambda x:
            x > 0.06*x.data.max()) 
    logger.info(strm("before filter",refl_blip_ranges))
    refl_blip_ranges = filter_range(refl_blip_ranges)  # repeats the filter range but for the reflected signal                
    refl_blip_ranges.sort(axis=0) # they are sorted by range size, not first/last
    logger.info(strm("after filter",refl_blip_ranges))
    assert refl_blip_ranges.shape[0] == 2, "seems to be more than two tuning blips "
    for thisrange in refl_blip_ranges:
        fl.plot(abs(d['ch',1]['t':tuple(thisrange)]), alpha=0.1, color='k',
                linewidth=10)
    # }}}
    
    # {{{ apply a linear phase to find any remaining fine offset of the pulse,
    #     and demodulate
    f_shift = nddata(r_[-0.1e6:0.1e6:200j],'f_test')
    # perform and store 200 frequency de-modulations of signal
    test_array = d['ch',0] * exp(-1j*2*pi*f_shift*d.fromaxis('t'))
    # performs frequency shift to control signal
    fine_adj_frq = test_array.sum('t').run(abs).argmax('f_test').item()
    # when modulating by same frequency of the waveform,
    # abs(sum(waveform)) will be a maximum
    center_frq += fine_adj_frq
    logger.info(strm(("found center frequency at %0.5f MHz"%(center_frq/1e6))))
    d.ft('t') #Fourier Transform into freq domain
    d.setaxis('t', lambda x: x - fine_adj_frq) #apply shift to x axis
    fl.next('frequency domain\n%s'%searchstr)
    fl.plot(d['t':(None,40e6)], label='demod and sliced',
            alpha=0.5)
    # }}}
    
    # {{{ use the "standard cost function" to determine the
    #     t=0 (treat decay as an FID)
    d.ift('t')
    first_blip = d['ch',1][
            't':tuple(refl_blip_ranges[0]+r_[-1e-6,1e-6])].C
    if standard_cost:
        first_blip.ft('t')
        fl.next('test time axis')
        t_shift = nddata(r_[-0.2e-6:0.2e-6:1000j]+pulse_range[0],
                't_shift').set_units('t_shift','s')
        # perform and store 1000 time shifts 
        test_data = first_blip * exp(
                -1j*2*pi*t_shift*first_blip.fromaxis('t'))
        test_data_ph = test_data.C.sum('t')
        test_data_ph /= abs(test_data_ph)
        test_data /= test_data_ph
        test_data.run(real).run(abs).sum('t')
        fl.plot(test_data,'.')
        # determine time zero
        time_zero = test_data.argmin('t_shift').item()
    else:
        time_zero = abs(first_blip).argmax('t').item()
    d.setaxis('t', lambda x: x-time_zero).register_axis({'t':0})
    refl_blip_ranges -= time_zero
    pulse_range -= time_zero
    d = d['t':(-10e6,10e6)] # slice out frequencies with signal
    #}}}
    
    #{{{zeroth order phase correction
    for j in range(2):
        fl.basename = "channel %d"%(j+1)
        ph0 = zeroth_order_ph(d['ch',j], fl=fl)
        d['ch',j] /= ph0
    fl.basename = None
    #}}}

    # {{{ 
    fl.next('after all corrections are complete')
    for j in range(2):
        fl.plot(d['ch',j].real,
                label='ch %d real'%(j+1), alpha=0.5)
        fl.plot(d['ch',j].imag,
                label='ch %d imag'%(j+1), alpha=0.5)
        fl.plot(abs(d['ch',j]), linewidth=3, color='k',
                label='ch %d abs'%(j+1), alpha=0.3)
    # }}}

    #{{{ to plot the transfer function, we need to pick an impulse
        # of finite width, or else we get a bunch of noise
    if show_transfer_func:
        transf_range = (-0.5e-6,3e-6)
        fl.next('the transfer function')
        impulse = exp(-d.fromaxis('t')**2/2/(0.03e-6)**2) #impulse function
        ## the following gives a possibility for a causal impulse
        fl.plot(impulse['t':transf_range], alpha=0.5, color='k', label='impulse')
        #plots impulse function in range of transfer function
        d.ft('t') #Fourier Transforms into freq domain
        transf = d['ch',1]/d['ch',0] #defining transfer function
        impulse.ft('t') #applies FT to impulse function
        response = impulse*transf #defines response
        response.ift('t') #Inverse Fourier transforms the response (which includes the impulse)
        response = response['t':transf_range] #defines x axis range of response
        fl.plot(response.real, alpha=0.5, label='response, real')
        fl.plot(response.imag, alpha=0.5, label='response, imag')
        fl.plot(abs(response), alpha=0.3, linewidth=3, label='response, abs')
    #}}}

    #{{{ fits curve to find Q
    dw = diff(d.getaxis('t')[0:2]).item()
    for thislabel,decay in [('initial',d['ch',1]['t':(refl_blip_ranges[0,0]-1e-6,refl_blip_ranges[1,0]-1e-6)]),
            ('final',d['ch',1]['t':(refl_blip_ranges[1,0]-1e-6,None)])]:
        max_t = abs(decay).argmax('t').item()
        decay = decay['t':(max_t,None)]
        decay = decay.setaxis('t',lambda x: x-decay.getaxis('t')[0])
        decay = decay['t':(None,2e6)] # none seem to exceed this -- just for plotting, etc
        fl.next('Plotting the decay slice for the %s blip'%thislabel)
        fl.plot(abs(decay), linewidth=3, alpha=0.3, color='k', label='starts at %g s'%max_t)
        fitfunc = lambda p: p[0]*exp(-decay.fromaxis('t')*p[1])+p[2] 
        #defines fit function as p0exp(-(t-t0)*p1)+p2
        p_ini = r_[decay['t',0].data.real.max(), 1/0.5e-6, 0] #why is there a third number (0) here?
        fl.plot(fitfunc(p_ini), ':', label='initial guess', alpha=0.5) 
        #applies the fit function to the initial point of the decay
        residual = lambda p: fitfunc(p).data.real - decay.data.real
        #subtracts the difference from the fit and the real data
        p_opt, success = leastsq(residual, p_ini[:])
        #fitting the data with least squares
        assert success > 0 & success < 5, "fit not successful"
        Q = 1./p_opt[1]*2*pi*center_frq 
        #relating the fit function to Q
        fl.plot(fitfunc(p_opt), label='fit, Q=%0.1f'%Q, alpha=0.5)
        fl.plot(decay.real, label='data (real)', alpha=0.5)
        fl.plot(decay.imag, label='data (imag, not fit)', alpha=0.5)
        #}}}
fl.show()
quit()

