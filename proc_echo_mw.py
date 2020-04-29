from pyspecdata import *
from Utility import dBm2power
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import * 
from sympy import symbols
rcParams["savefig.transparent"] = True
fl = fl_mod()
t2 = symbols('t2')
# leave this as a loop, so you can load multiple files
for searchstr,freq_range,time_range in [
        ["200306_DNP_lg_probe_w34.*",(-300,300),(None,0.05)]
        ]:

    s = load_data(searchstr,'ODNP',None)

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
