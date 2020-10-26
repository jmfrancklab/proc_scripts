from pyspecdata import *
from scipy.optimize import minimize,leastsq
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
from scipy.signal import tukey
#do_slice should be true when using AFG, and false when using spincore
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
        #('200103_pulse_1','test_equip','capture1','square_wave_capture_v1',True),
        #('201020_sol_probe_1','test_equip','capture1','square_wave_capture_v1',True)
        #('201009_coilE_4','test_equip','capture1','square_wave_capture_v1',True),
        ('201026_sqwv_cap_probe_1','test_equip','capture1','square_wave_capture_v1',True)
        ]:
    d = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict) 
    print(d.getaxis('t'))
    fl.next('c')#'Raw signal %s'%searchstr)
    fl.plot(d['ch',0], alpha=0.5, label='control') # turning off human units forces plot in just V
    fl.plot(d['ch',1], alpha=0.5, label='reflection')
    #fl.show();quit()
    # }}}
    
    # {{{ determining center frequency and convert to
    # analytic signal, show analytic signal
    d.ft('t')#,shift=True) #Fourier Transform into freq domain
    d = 2*d['t':(0,None)] # throw out negative frequencies and low-pass (analytic signal -- 2 to preserve magnitude)
    #to negated the "1/2" in "(1/2)(aexp[iwt]+a*exp[-iwt])
    # ALSO -- everything should be less than 40 MHz for sure
    if do_slice:
        d = d['t':(0,50e6)]
        d[lambda x: abs(x) < 1e-10] = 0
    else:
        d['t':(0,3e6)] = 0 # filter out low-frq noise, which becomes high-frq noise on demod
        d[lambda x: abs(x) < 1e-10] = 0
    center_frq = abs(d['ch',0]).argmax('t').item() # the center frequency is now the max of the freq peak    
    logger.info(strm(("initial guess at center frequency at %0.5f MHz"%(center_frq/1e6))))
    logger.info(strm(center_frq))
    freq_guess = abs(d['ch',1]).argmax('t').item()
    freq_range = r_[-10,10]*1e6 +freq_guess
    d['t':(0,freq_range[0])] = 0
    d['t':(freq_range[1],None)] = 0
    tukey_filter = d.fromaxis('t')['t':tuple(freq_range)].run(lambda x: tukey(len(x)))
    d['t':tuple(freq_range)] *= tukey_filter
    fl.next('showing freq dist.')
    for j in d.getaxis('ch'):
        #fl.plot(abs(forplot)['ch':j],alpha=0.5,plottype='semilogy',label=f'CH{j} orig')
        fl.plot(abs(d)['ch':j][lambda x: abs(x) > 1e-10],alpha=0.5,plottype='semilogy',label=f'CH{j} filtered')
        fl.grid()
    fl.next('frequency domain\n%s'%searchstr)
    fl.plot(abs(d['t':(None,40e6)]),label='Raw signal in freq domain,\nshows a bandwidth of about 20 MHz', alpha=0.5)
    axvline(x=center_frq/1e6)
    d.ift('t') #Inverse Fourier Transform back to time domain to display the decaying exponential
        # {{{ determine the frequency from the phase gradient during the pulse
    dt = diff(d.getaxis("t")[r_[0, 1]]).item()
    pulse_slice = d["ch", 0].contiguous(lambda x: abs(x) > 0.5*abs(x).data.max())[0] #defines pulse slice based on control signal
    d.setaxis("t", lambda x: x - pulse_slice[0]).register_axis({"t": 0}) #resets t axis around pulse slice
    pulse_slice -= pulse_slice[0]
    #{{{ Not sure what this portion is doing...is this similar to a time shift? or first order phase correction?
    d = d["t" : tuple(pulse_slice + r_[-0.5e-6, 2e-6])]
    pulse_middle = d["ch", 0]["t" : tuple(pulse_slice + r_[+0.5e-6, -0.5e-6])]
    ph_diff = pulse_middle["t", 1:] / pulse_middle["t", :-1]
    ph_diff.sum("t")
    ph_diff = ph_diff.angle.item()
    frq = ph_diff/dt/2/pi
    #}}}
    # }}}
    print("frq:", frq)
    d *= exp(-1j*2*pi*frq*d.fromaxis("t")) #convolution
    ph0 = d["ch", 0].C.sum("t").item() #pseudo 0th order phase correction
    ph0 /= abs(ph0)
    d /= ph0
    fl.next('Absolute value of analytic signal, %s'%searchstr)
    fl.plot(abs(d['ch',0]), alpha=0.5, label='control') #plot the 'envelope' of the control 
    fl.plot(abs(d['ch',1]), alpha=0.5, label='reflection') #plot the 'envelope' of the reflection so no more oscillating signal
    #fl.show();quit()
    #{{{determine the start and stop points for both the pulse, as well as the two tuning blips
    scalar_refl = d["ch", 1]["t":(2e-6, 8e-6)].mean("t").item()
    blip_range = r_[-0.13e-6, 1.5e-6] #defining decay slice
    first_blip = -d["ch", 1:2]["t" : tuple(blip_range)] + scalar_refl #correcting first blip
    #{{{ doing 0th order correction type thing again? why? we did this in lines 99-101...
    ph0_blip = first_blip["t", abs(first_blip).argmax("t", raw_index=True).item()]
    ph0_blip /= abs(ph0_blip)
    d1 = first_blip/ph0_blip
    fl.next('smoothed decay of first and second blip')
    fl.plot(d1.real,label='first blip real')
    fl.plot(d1.imag,label='first blip imag')
    fl.plot(abs(d1),label='first blip abs')
    fl.twinx(orig=False)
    ax = gca()
    gridandtick(ax)
    ax.grid(False)
    fl.twinx(orig=True)
    fl.grid
    second_blip = d['ch', 1:2]['t':tuple(blip_range + pulse_slice[1])].setaxis('t',lambda x: x - pulse_slice[1])
    d2 = second_blip/ph0_blip
    fl.plot(d2.real,label='second blip real')
    fl.plot(d2.imag,label='second blip imag')
    fl.plot(abs(d2),label='second blip abs')
    #fl.show();quit()
    #{{{ fits curve to find Q
    print(ndshape(d1))
    print(ndshape(d2))
    decay = (abs(d1)+abs(d2))/2
    fl.next('decay')
    fl.plot(decay)
    #decay_start = decay.argmax('t').item()
    #decay = decay['t':(decay_start,None)]
    decay = decay['t':(69e-9,1200)]
    fl.next('Plotting the decay slice')
    fl.plot(decay, linewidth=3, alpha=0.3, color='k')
    print(decay.getaxis('ch'))
    decay = decay['ch',0]
    print(ndshape(decay))
    f = fitdata(decay)
    A,B,C,t = symbols("A B C t",real=True)
    f.functional_form = A*e**(-t*B)
    fl.next('fit')
    fl.plot(decay,'o',label='data')
    f.fit()
    f.set_units('t','ns')
    print("output:",f.output())
    print("latex:",f.latex())
    Q = 1./f.output('B')*2*pi*center_frq
    fl.plot(f.eval(100).set_units('t','s'),label='fit, Q=%0.1f'%Q)
fl.show()
quit()

