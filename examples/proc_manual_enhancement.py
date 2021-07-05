from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from sympy import symbols, Symbol, latex,limit,init_printing
from numpy import *
import matplotlib.pyplot as plt
#{{{input parameters
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
fl = fl_mod()
measured_vs_actual = 22. # how many dB down the split + measured power is from
#                          the forward power
#}}}

s = find_file('210702_100mM_TEMPO_hexane_test_2',exp_type='ODNP_NMR_comp/test_equipment', expno='enhancement_curve',
        postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=fl)
print(s.get_prop('acq_params'))
fl.next('raw data')
fl.image(s)
#fl.show();quit()
s = s['t2':(-5e3,5e3)]
ph0 = s['power',-4].sum('t2')
ph0 /= abs(ph0)
s /= ph0
fl.next('phased')
fl.image(s)
#fl.show();quit()
s = s['t2':(-300,-55)]
fl.next('sliced to integration bounds')
fl.image(s)
s = s['ph1',1]
fl.next('sliced to integration bounds\nselect coh. pathway\nreal')
fl.image(s.real)
s.integrate('t2')
s /= max(s.data.real)
#{getting power axis
print(s.get_prop('acq_params')['meter_powers'])
s.setaxis('power',r_[-9999,array(s.get_prop('acq_params')['meter_powers'])])
print("here are the dBm",s.getaxis('power'))
s.setaxis('power', lambda x:
        1e-3*10**((x+measured_vs_actual)/10.))
print("here are the powers",s.getaxis('power'))
s.set_units('power','W')
#}
fl.next('simple integral')
fl.plot(s['power',:-3],'ko', human_units=False) # human_units = False because the range for the red points tries to force mW, which is incompatible with W
fl.plot(s['power',-3:],'ro', human_units=False)
plt.ylabel('integrated $^1$H NMR signal')
#{Enhancement
s /= s['power',0]
fl.next('enhancement curve')
fl.plot(s['power',:-3],'ko', human_units=False)
fl.plot(s['power',-3:],'ro', human_units=False)
fl.show()
#}
