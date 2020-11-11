from pyspecdata import *
from scipy.optimize import minimize,leastsq
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
from scipy.signal import tukey
do_slice = False # slice frequencies and downsample -- in my hands, seems to decrease the quality of the fit 
standard_cost = False # use the abs real to determine t=0 for the blip -- this actually doesn't seem to work, so just use the max
show_transfer_func = False # show the transfer function -- will be especially useful for processing non-square shapes
logger = init_logging('info')
#init_logging(level='debug')
# 2to3 JF 1/31

fl = figlist_var()
 # {{{ load data, set units, show raw dataz
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201106_chirp_shortingcap_1','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    a = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=postproc_dict)
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201030_chirp_RM_probe_1','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    b = find_file(searchstr,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict)
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201106_chirp_50ohm_1','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    c = find_file(searchstr,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict)
for searchstr,exp_type,nodename,postproc,corrected_volt in [
        ('201106_chirp_empty_1','ODNP_NMR_comp','capture1','chirp',True)
        ]:
    d = find_file(searchstr,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict)
fl.next('Raw signal')
fl.plot(a['ch',0], alpha=0.5, label='shorting cap control')    
fl.plot(a['ch',1], alpha=0.5, label='shorting cap refl')
fl.plot(b['ch',0], alpha=0.5, label='RM probe control')    
fl.plot(b['ch',1], alpha=0.5, label='RM probe refl')
fl.plot(c['ch',0], alpha=0.5, label='50 ohm control')    
fl.plot(c['ch',1], alpha=0.5, label='50 ohm refl')
fl.plot(d['ch',0], alpha=0.5, label='open control')    
fl.plot(d['ch',1], alpha=0.5, label='open refl')

a.ft('t',shift=True)
b.ft('t',shift=True)
c.ft('t',shift=True)
d.ft('t',shift=True)
fl.next('freq domain for RM probe with control')
fl.plot(a['ch',0],alpha=0.5,label='control shorting cap')
fl.plot(a['ch',1], alpha=0.5, label='reflection shorting cap')
fl.plot(b['ch',0],alpha=0.5,label='control RM probe')
fl.plot(b['ch',1],alpha=0.5,label='reflection RM probe')
fl.plot(c['ch',0],alpha=0.5,label='control 50 ohm')
fl.plot(c['ch',1],alpha=0.5,label='reflection 50 ohm')
fl.plot(d['ch',0],alpha=0.5,label='control open')
fl.plot(d['ch',1],alpha=0.5,label='reflection open')
for j in range(ndshape(a)['ch']):
    a = a['t':(10e6,18e6)]
a_control = a['ch',0]['t':(0,None)]
a_refl = a['ch',1]['t':(0,None)]
a = a_refl/a_control
for j in range(ndshape(b)['ch']):
    b = b['t':(10e6,18e6)]
b_control = b['ch',0]['t':(0,None)]
b_refl = b['ch',1]['t':(0,None)]
b = b_refl/b_control
for j in range(ndshape(c)['ch']):
    c = c['t':(10e6,18e6)]
c_control = c['ch',0]['t':(0,None)]
c_refl = c['ch',1]['t':(0,None)]
c = c_refl/c_control
for j in range(ndshape(d)['ch']):
    d = d['t':(10e6,18e6)]
d_control = d['ch',0]['t':(0,None)]
d_refl = d['ch',1]['t':(0,None)]
d = d_refl/d_control

fl.next('when divided by control for %s'%searchstr)
fl.plot(abs(a),label='shorting cap')
fl.plot(abs(b),label='RM probe')
fl.plot(abs(c),label='50 ohm')
fl.plot(abs(d),label='open')
ylim(0,1.25)
fl.next('phase for controls')
fl.plot(a.angle,label='shorting cap')
fl.plot(b.angle,label='RM probe')
fl.plot(c.angle,label='50 ohm')
fl.plot(d.angle,label='open')
fl.show();quit()

