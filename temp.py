from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from sympy import symbols, Symbol, latex, limit, init_printing
from numpy import *
plt.rcParams.update({
    "figure.facecolor": (1.,1.,1.,0.), # clear
    "axes.facecolor": (1.,1.,1.,0.9), # 90% transparent white
    "savefig.facecolor": (1.,1.,1.,0.), #clear
    })
logger = init_logging("info")
t2 = symbols('t2')
fl = fl_mod()

s = find_file('210617_T177R1a_pR_DDM_ODNP', exp_type='odnp',expno='FIR_0dBm',
        postproc='spincore_IR_v1', lookup=postproc_dict, fl=None)
s.reorder(['ph1','ph2','vd','t2'])
fl.next('raw data')
fl.image(s)
fl.next('raw data - time domain')
fl.image(s.C.ift('t2')['t2':(0,100e-3)])
s = s['t2':(-1e3,1e3)]
fl.next('sliced to integration bounds')
fl.image(s)
#{phasing
s.ift('t2')
s.ift(['ph1','ph2'])
phasing = s['t2',0].C
phasing.data *= 0
phasing.ft(['ph1','ph2'])
phasing['ph1',0]['ph2',1] = 1
phasing.ift(['ph1','ph2'])
ph0 = s['t2':0]/phasing
ph0 /= abs(ph0)
s /= ph0
s.ft(['ph1','ph2'])
s.ft('t2')
fl.next('phased')
fl.image(s)
#}
d = s['ph1',1].C
d.integrate('t2')
d /= max(d.data)
fl.next('simple integration')
fl.plot(d)
#{correlation alignment
#}
fl.show()
