from pyspecdata import *
from Utility import dBm2power
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import * 
from sympy import symbols
rcParams["savefig.transparent"] = True
t2 = symbols('t2')

s = load_data("200306_DNP_lg_probe_w34.*")
freq_range = (-300,300)
time_range = (None,0.05)

prog_power = s.getaxis('power').copy()
print("programmed powers",prog_power)
s.setaxis('power',r_[
    0,dBm2power(array(s.get_prop('meter_powers'))+20)]
    ).set_units('power','W')
print("meter powers",s.get_prop('meter_powers'))
print("actual powers",s.getaxis('power'))
print("ratio of actual to programmed power",
           s.getaxis('power')/prog_power)
nPoints = s.get_prop('acq_params')['nPoints']
SW_kHz = s.get_prop('acq_params')['SW_kHz']
nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
s.set_units('t','s')
s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
s.labels({'ph2':r_[0.,2.]/4,
    'ph1':r_[0.,1.,2.,3.]/4})
s.reorder(['ph2','ph1'])
s.ft('t2',shift=True)
s.ft(['ph1','ph2'])
fl.next('all data: frequency domain')
fl.image(s)
fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
        s,freq_range)
s = s['t2':freq_range]
s.ift('t2')
residual,best_shift = hermitian_function_test(s[
    'ph2',-2]['ph1',1])
fl.next('hermitian test')
fl.plot(residual)
s.setaxis('t2',lambda x: x-best_shift)
s.register_axis({'t2':0}, nearest=False)
    # {{{ implement zeroth-order correction
    # note that it's not going to have only one
    # dimension, b/c we will have at least a power
    # dimension

s = FID(s)
#ph0 = s['t2':0]['ph2',-2]['ph1',1]
#s /= zeroth_order_ph(ph0, fl=fl)
#if s['t2':0]['ph2',-2]['ph1',1]['power',0].real < 0:
#   s *= -1
# }}}
#fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range',
#        s,time_range)
#s = s['t2':time_range]
    # {{{ all of this is to check and see if we think
    # we can add the two sides of the echo to increase
    # our SNR -- right now, it doesn't look like it
#fl.next('echo mirror test')
#echo_start = s.getaxis('t2')[0]
#dw = diff(s.getaxis('t2')[r_[0,1]]).item()
#centered_echo = s['t2':(echo_start,-echo_start+dw)]
#plotdata = abs(centered_echo/centered_echo['t2',::-1])
#plotdata[lambda x: x>2] = 2
#fl.image(plotdata)
## }}}
#fl.next('apodize and zero fill')
#R = 5.0/(time_range[-1]) # assume full decay by end time
#s *= exp(-s.fromaxis('t2')*R)
#s.ft('t2',pad=1024)
#fl.image(s)
#   # {{{ select FID
#s.ift('t2')
#s = s['ph2',-2]['ph1',1]['t2':(0,None)]
#s['t2',0] *= 0.5
#s.ft('t2')
## }}}
fl.next('compare highest power to no power')
idx_maxpower = argmax(s.getaxis('power'))
fl.plot(s['power',0])
fl.plot(s['power',idx_maxpower])
fl.next('full enhancement curve')
fl.plot(s)
fl.next('enhancement')
enhancement = s['t2':(-50,50)].C.sum('t2').real
enhancement /= enhancement['power',0]
fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)
ylabel('Enhancement')
fl.show();quit()
