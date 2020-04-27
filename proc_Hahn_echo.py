from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import load_data
fl = figlist_var()
s = load_data('191031_echo_4_5.h5')
nEchoes = 1
nPhaseSteps = 1
SW_kHz = 24.0
nPoints = 1024
s.set_units('t','s')
print(s.get_prop('acq_params'))
print(s.get_prop('nScans'))
fl.next('raw data')
fl.plot(s.real,alpha=0.4)
#fl.plot(s.imag,alpha=0.4)
#fl.plot(abs(s),':',c='k',alpha=0.4)
s.ft('t',shift=True)
fl.next('comp raw data - FT')
fl.plot(s.real,alpha=0.4)
fl.plot(s.imag,alpha=0.4)
#fl.plot(abs(s),c='red')
#fl.next('comp raw data - FT')
#fl.plot(s.real,alpha=0.4)
#fl.plot(s.imag,alpha=0.4)
#fl.plot(abs(s),':',c='k',alpha=0.4)
fl.show()
