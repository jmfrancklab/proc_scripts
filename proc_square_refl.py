from pyspecdata import *
from scipy.optimize import minimize,leastsq

from hermitian_function_test import zeroth_order_ph

init_logging(level='debug')
# 2to3 JF 1/31

fl = figlist_var()
 # {{{ load data, set units, show raw data
for date, id_string,corrected_volt in [
        #('180806','pulse_reflection',True),
        #('181001','sprobe_t2',True),
        #('181001','sprobe_t4',True),
        #('181103','probe',True),
        ('200110','pulse_2',True),
        #('200312','chirp_coile_4',True),
        #('200103','pulse_1',True),
        ]:
    d = nddata_hdf5(date+'_'+id_string+'.h5/capture1',
                directory=getDATADIR(exp_type='test_equip'))
    d.set_units('t','s').name('Amplitude').set_units('V')
    print(ndshape(d))
    fl.next('Raw signal %s'%id_string)
    if date == '191212':
        d['ch',0] *= 0.5
    fl.plot(d['ch',0], alpha=0.5, label='control') # turning off human units forces plot in just V
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    # }}}
    # {{{ determining center frequency and convert to
    # analytic signal, show analytic signal
    d.ft('t',shift=True) #Fourier Transform into freq domain
    d = 2*d['t':(0,None)] # throw out negative frequencies and low-pass (analytic signal -- 2 to preserve magnitude)
    #to negated the "1/2" in "(1/2)(aexp[iwt]+a*exp[-iwt])
    center_frq = abs(d['ch',0]).argmax('t').item() # the center frequency is now the max of the freq peak    
    print(("initial guess at center frequency at %0.5f MHz"%(center_frq/1e6)))
    print(center_frq)
    fl.next('Raw signal in freq domain, abs: %s'%id_string)
    fl.plot(abs(d), human_units=False)
    axvline(x=center_frq)
    d.ift('t') #Inverse Fourier Transform back to time domain to display the decaying exponential
    fl.next('Absolute value of analytic signal, %s'%id_string)
    fl.plot(abs(d['ch',0]), alpha=0.5, label='control') #plot the 'envelope' of the control 
    fl.plot(abs(d['ch',1]), alpha=0.5, label='reflection') #plot the 'envelope' of the reflection so no more oscillating signal
    # }}}
    # {{{ determine the start and stop points for both
    # the pulse, as well as the two tuning blips
    pulse_range = abs(d['ch',0]).contiguous(lambda x:  # returns list of limits for which the lambda function holds true
            x > 0.5*x.data.max())                      # So will define pulse_range as all x values where the signal goes 
                                                       # above half the max of the signal
    #pulse_range = array(j for j in pulse_range if
    #        diff(j).item() > 0.1e-6)
    # RATHER -- use numpy
    def filter_range(thisrange):                      
        mask = diff(thisrange, axis=1) > 0.1e-6 * ones((1,2))  #filters out when the signal only goes 0.1-0.2 us above half max
        thisrange = thisrange[mask].reshape((-1,2))
        return thisrange
    pulse_range = filter_range(pulse_range)
    if not pulse_range.shape[0] == 1:
        print(("seems to be more than one pulse -- on starting at " 
                + ','.join(('start '+str(j[0])+' length '+str(diff(j)) for j in pulse_range))))   # If there is more than one section that goes above half max
# it assumes theres more than one pulse 
      
    pulse_range = pulse_range[0,:]
    fl.plot(abs(d['ch',0]['t':tuple(pulse_range)]), alpha=0.1, color='k',  #shades in the section of pulse range (above half max) for 
            linewidth=10)                                                  #control 
    refl_blip_ranges = abs(d['ch',1]).contiguous(lambda x:
            x > 0.05*x.data.max()) 
    refl_blip_ranges = filter_range(refl_blip_ranges)                      # repeats the filter range but for the reflected signal                
    assert refl_blip_ranges.shape[0] == 2, "seems to be more than two tuning blips "
    for thisrange in refl_blip_ranges:
        fl.plot(abs(d['ch',1]['t':tuple(thisrange)]), alpha=0.1, color='k',
                linewidth=10)
        #fl.show();quit()
    # }}}
    # {{{ apply a linear phase to find the frequency of the pulse (control)
    f_shift = nddata(r_[-0.1e6:0.1e6:200j]+center_frq,'f_test')
    # perform and store 200 frequency de-modulations of signal
    test_array = d['ch',0] * exp(-1j*2*pi*f_shift*d.fromaxis('t'))
    # performs frequency shift to control signal
    center_frq = test_array.sum('t').run(abs).argmax('f_test').item()
    # when modulating by same frequency of the waveform,
    # abs(sum(waveform)) will be a maximum
    print(("found center frequency at %0.5f MHz"%(center_frq/1e6)))
    # }}}
    # {{{ creates demodulated reflection 
    d.ft('t') #Fourier Transform into freq domain
    d.setaxis('t', lambda x: x - center_frq) #apply shift to x axis
    fl.next('test time axis')
    t_shift_testvals = r_[-1e-6:1e-6:1000j]+pulse_range[0]
    # perform and store 1000 frequency de-modulations of signal
    t_shift = nddata(t_shift_testvals,'t_shift')
    # perform and store 1000 time shifts 
    test_data = d['ch',0] * exp(-1j*2*pi*t_shift*d.fromaxis('t'))
    #apply t shift to control data
    test_data_ph = test_data.C.sum('t')
    test_data_ph /= abs(test_data_ph)
    test_data /= test_data_ph
    test_data.run(real).run(abs).sum('t')
    fl.plot(test_data,'.')
    #fl.show();quit()
    #Cost function to correct zeroth order phasing
    pulse_start = test_data.argmin('t_shift').data
    #all data is now starting at min of t
    # TODO here -- filter out frequencies very far from the center frequency
    d.ift('t') #Inverse Fourier Transform into t domain 
    # }}}
    for j in range(2):
        ph0 = zeroth_order_ph(d['ch',j])
        d['ch',j] /= ph0
    fl.next('demodulated and phased data')
    #d = d['t':(2e-06,4e-06)]

    for j in range(2):
        fl.plot(d['ch',j], label='channel %d real'%j,
                alpha=0.5)
        #fl.plot(d['ch',j].imag, label='channel %d imag'%j,
                #alpha=0.5)
        #fl.show();quit()
    d.setaxis('t',lambda x: x-pulse_start)
    #fl.show();quit()
    #print("NOTE!!! the demodulated reflection looks bad -- fix it")
    # to use the phase of the reference to set both, we could do:
    # pulse_phase = d['ch',0].C.sum('t')
    # but I don't know if that's reasonable -- rather I just phase both independently:
    #pulse_phase = d.C.sum('t')
    #pulse_phase /= abs(pulse_phase)
    #d /= pulse_phase

    #cost function for phase correction
    for j,l in enumerate(['control','reflection']):
        fl.next('adjusted analytic '+l)
        fl.plot(d['ch',j].real, alpha=0.3, label='real')
        fl.plot(d['ch',j].imag, alpha=0.3, label='imag')
        fl.plot(abs(d['ch',j]), alpha=0.3, color='k', linewidth=2,
                label='abs')
    # }}}
    # {{{ to plot the transfer function, we need to pick an impulse
    # of finite width, or else we get a bunch of noise
    transf_range = (-0.5e-6,3e-6)
    fl.next('the transfer function')
    impulse = exp(-d.fromaxis('t')**2/2/(0.03e-6)**2) #impulse function
    ## the following gives a possibility for a causal impulse
    #impulse = exp(-abs(d.fromaxis('t'))/0.01e-6)
    #impulse['t':(None,0)] = 0
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
    #fl.show();quit()
    fl.next('Plotting the decay slice')

    d.ift('t')

    #print(nddata(refl_blip_ranges))

    dw = diff(d.getaxis('t')[0:2]).item()
    def phasecorrect(d):
        fl.push_marker() 
        ph1 = nddata(r_[-5:5:70j]*dw,'phcorr')
        dx = diff(ph1.getaxis('phcorr')[r_[0,1]]).item()
        ph1 = exp(-1j*2*pi*ph1*d.fromaxis('t'))
        d_cost = d * ph1
        ph0 = d_cost.C.sum('t')
        ph0 /= abs(ph0)
        d_cost /= ph0
        fl.next('phasing cost function')
        d_cost.run(real).run(abs).sum('t')
        fl.plot(d_cost,'.')
        ph1_opt = d_cost.argmin('phcorr').item()
        print('optimal phase correction',repr(ph1_opt))
        # }}}
        # {{{ apply the phase corrections
        def applyphase(arg,ph1):
            arg *= exp(-1j*2*pi*ph1*arg.fromaxis('t'))
            ph0 = arg.C.sum('t')
            ph0 /= abs(ph0)
            arg /= ph0
            return arg
        def costfun(ph1):
            if type(ph1) is ndarray:
                ph1 = ph1.item()
            temp = d.C
            retval = applyphase(temp,ph1).run(real).run(abs).sum('t').item()
            return retval
        print("rough opt cost function is",costfun(ph1_opt))
        r = minimize(costfun,ph1_opt,
                bounds=[(ph1_opt-dx,ph1_opt+dx)])
        assert r.success
        d = applyphase(d,r.x.item())
        fl.plot(r.x,r.fun,'x')
        fl.pop_marker()
        return d
    d.ft('t')
    d = phasecorrect(d['ch',1])
    fl.next('phased f domain')
    fl.plot(d, label = 'refl real')
    #fl.show();quit()
        #d.ft('t')
    d['t':(30,None)]=0
    #d['t':(None,-20)]=0
    fl.next('phased t domain')
    d.ift('t')
    fl.plot(d)
    #fl.show();quit()
    #fl.plot(d['ch',1].imag, label = 'refl imag')
    #fl.plot(d['ch',0].real, label = 'control real')
    #fl.plot(d['ch',0].imag, label = 'control imag')
    #fl.show();quit()
    #quit()# Inverse Fourier Transform into t domain
    fl.next('decay')

    decay = d['t':(0.5e-6,2e-6)]
        #0.5*(refl_blip_ranges[0,1]+refl_blip_ranges[1,0]))]
    fl.plot(decay)
    #fl.show();quit()
    # slice out a range from the start of the first
    # blip up to halfway between the END of the first
    # blip and the start of the second
    fl.next('fit')
    decay = decay.setaxis('t',lambda x: x-decay.getaxis('t')[0])
    # }}}
    fitfunc = lambda p: p[0]*exp(-decay.fromaxis('t')*p[1])+p[2] 
    #defines fit function as p0exp(-(t-t0)*p1)+p2
    p_ini = r_[decay['t',0].data.real.max(),1/0.5e-6,0] #why is there a third number (0) here?
    fl.plot(fitfunc(p_ini), ':', label='initial guess') 
    #applies the fit function to the initial point of the decay
    residual = lambda p: fitfunc(p).data.real - decay.data.real
    #subtracts the difference from the fit and the real data
    p_opt, success = leastsq(residual, p_ini[:])
    #fitting the data with least squares
    assert success > 0 & success < 5, "fit not successful"
    Q = 1./p_opt[1]*2*pi*center_frq 
    #relating the fit function to Q
    fl.plot(fitfunc(p_opt), label='fit')
    fl.plot(decay, label='data')
    fl.plot(decay.imag, label='data (imag, not fit)')
    print(Q)
fl.show()
