from pyspecdata import *
from Utility import dBm2power
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import * 
from sympy import symbols
rcParams["savefig.transparent"] = True
fl = fl_mod()
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
s = FID(s,(None,0.05))
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
