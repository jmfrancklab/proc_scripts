#{Imports & Initialization
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from sympy import symbols, Symbol, latex,limit,init_printing
from numpy import *
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.0),  # clear ## 1,1,1,0.9 == 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
fl = fl_mod()
measured_vs_actual = 22. # how many dB down the split + measured power is from
                         # the forward power
#}
#{Input parameters & Load Data
f_range =(-250,75)         #(-250,100)#(-200,200) # Frequency cutoff in Hz. 
outname = '210714_A174R1a_pR_DHPC_enhancement'#DDM_enhancement'#'210714_150uM_TEMPOL_enhancement'
# If you don't want a csv of enhancement data, put None.
s = find_file('210714_A174R1a_pR_DHPC_ODNP',exp_type='odnp',expno='enhancement',#'enhancement1',
        postproc='spincore_ODNP_v2',lookup=postproc_dict,fl=fl)
#s = find_file('210616_S175R1a_pR_DDM_ODNP',exp_type='odnp',expno='enhancement',
#        postproc='spincore_ODNP_v2',lookup=postproc_dict,fl=fl)
#s = find_file('210617_T177R1a_pR_DDM_ODNP',exp_type='odnp', expno='enhancement',
#    postproc='spincore_ODNP_v2',lookup=postproc_dict,fl=fl)
#s = find_file('210507_TEMPOL_150uM__cap_probe_DNP_1',exp_type='odnp',expno='signal',
#        postproc='spincore_ODNP_v2',lookup=postproc_dict,fl=fl)
#}
fl.next('raw data')
s = s['t2':(-1e3,1e3)]
fl.image(s)
ph0 = s['power',-4].sum('t2')
ph0 /= abs(ph0)
s /= ph0
fl.next('phased')
fl.image(s)
s = s['t2':f_range]
fl.next('sliced to integration bounds')
fl.image(s)
if 'ph2' in s.dimlabels:
    s = s['ph1',1]['ph2',0]
else:
    s = s['ph1',1]
fl.next('sliced to integration bounds\nselect coh. pathway\nreal')
fl.image(s.real)
s.integrate('t2')
s /= max(s.data.real)
#{getting power axis
s.setaxis('power',r_[-9999,array(s.get_prop('meter_powers'))])
print("here are the dBm",s.getaxis('power'))
s.setaxis('power', lambda x:
        1e-3*10**((x+measured_vs_actual)/10.))
print("here are the powers",s.getaxis('power'))
s.set_units('power','W')
#}
#{Integral
fl.next('simple integral')
fl.plot(s['power',:-3],'ko', human_units=False) # human_units = False because the range for the red points tries to force mW, which is incompatible with W
fl.plot(s['power',-3:],'ro', human_units=False)
plt.ylabel('integrated $^1$H NMR signal')
#}
#{Enhancement
#s *= (207.4/113.1)
s /= s['power',0]
fl.next('enhancement curve')
fl.plot(s['power',:-3],'ko', human_units=False)
fl.plot(s['power',-3:],'ro', human_units=False)
fl.show()
#}
# { Exporting data from the enhancement curve as csv
if outname is not None:
    import pandas as pd
    reenhancementdf = pd.DataFrame(data=s.real.data,columns=['Re[E(p)]'])
    imenhancementdf = pd.DataFrame(data=s.imag.data,columns=['Im[E(p)]'])
    enhancementdf = pd.DataFrame(data=s.data,columns=['enhancement'])
    powerdf = pd.DataFrame(data=s.getaxis('power'),columns=['power'])
    enhancement = powerdf.join(reenhancementdf)
    enhancement2 = enhancement.join(imenhancementdf)
    enhancement3 = enhancement2.join(enhancementdf)
    enhancement3.to_csv('%s.csv'%outname,index=False)
#}
