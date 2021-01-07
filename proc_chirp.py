from pyspecdata import *
from scipy.optimize import minimize,leastsq
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
from scipy.signal import tukey
from pylab import *
do_slice = False # slice frequencies and downsample -- in my hands, seems to decrease the quality of the fit 
standard_cost = False # use the abs real to determine t=0 for the blip -- this actually doesn't seem to work, so just use the max
show_transfer_func = False # show the transfer function -- will be especially useful for processing non-square shapes
logger = init_logging('info')
#init_logging(level='debug')
# 2to3 JF 1/31

fl = figlist_var()
 # {{{ load data, set units, show raw dataz
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201113_triwave_control_2','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    a = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict)
fl.next('Raw signal')
fl.plot(a['ch',0], alpha=0.5, label='control')    
fl.plot(a['ch',1], alpha=0.5, label='control')
a.ft('t',shift=True)
fl.next('freq domain for RM probe with control')
fl.plot(a['ch',0],alpha=0.5,label='control')
fl.plot(a['ch',1], alpha=0.5, label='control')
a_control = a['ch',0]['t':(0,None)]
a_refl = a['ch',1]['t':(0,None)]
a = a_refl/a_control
#{{{idealized resonator
#intrinsic_Q = 100.
#f_0 = 14.8e6
#omega_0 = f_0*2*pi
#match_reactance = -50.0 # this is just based by eye off the plot, this is what I need to get the phase to cross zero at the characteristic impedance that I want
#LC = (omega_0)**-2
#L = sqrt(LC)
#C = sqrt(LC)
#C_match = -1.0/(match_reactance*omega_0)
#if C_match < 0:
#    raise ValueError("You chose your match capacitance to be negative -- you can't do this")
#delta_F = f_0/intrinsic_Q
#omega = r_[f_0-4*delta_F:f_0+4*delta_F:500j]*2*pi
## Q = R*sqrt(C/L) fo parallel, so
#R = intrinsic_Q/sqrt(C/L)
##parallel = 1./(1./(1j*omega*L) + 1j*omega*C)
#parallel_admittance = 1.0/(1j*omega*L) + (1j*omega*C) + 1./R
#parallel_impedance = 1.0/parallel_admittance
#tank_impedance = parallel_impedance + 1.0/(1j*omega*C_match)
#parallel_impedance = nddata(parallel_impedance,[-1],['f']).labels('f',omega/2/pi).set_units('f','Hz')
#tank_impedance = nddata(tank_impedance,[-1],['f']).labels('f',omega/2/pi).set_units('f','Hz')
#lorentzian_impedance = (omega_0/2)/(omega_0/intrinsic_Q/2+1j*(omega-omega_0)) # for reasons that I don't understand, this rate is in radians
#lorentzian_impedance = nddata(lorentzian_impedance,[-1],['f']).labels('f',omega/2/pi).set_units('f','Hz')
#
#impedance_at_match = 50.
#lorentzian_impedance /= intrinsic_Q
#lorentzian_impedance *= impedance_at_match

#Gamma_lorentzian = (lorentzian_impedance - 50.0)/(lorentzian_impedance + 50.0)
#Gamma_tank = (tank_impedance - 50.0)/(tank_impedance + 50.0)

fl.next('when divided by control for %s'%searchstr)
fl.plot(abs(a),label='shorting cap')
#fl.plot(abs(Gamma_lorentzian),'k',label = 'scaled lorentzian')
#fl.plot(abs(Gamma_tank),'k',alpha = 0.5,label = 'matched tank')

ylim(0,1.25)
#fl.next('phase for controls')
#fl.plot(a.angle,label='shorting cap')
#fl.plot(Gamma_lorentzian.angle,label='scaled lorentzian')
#fl.plot(Gamma_tank.angle,label='matched tank')
fl.show();quit()

