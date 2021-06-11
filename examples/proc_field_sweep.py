"""
Field Sweep
===========

A field sweep to detect the ESR via NMR-ODNP!
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from scipy.optimize import leastsq,minimize,basinhopping
from sympy import symbols
import pywt
fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 20e3
filename = '210611_S175R1a_pR_DDM_field_dep'
for nodename,postproc,label_str,freq_slice,field_slice in [
        ('32dBm_finer','field_sweep','Sams field sweep',(-700,700),(-400,300)),
        ]:
    s = find_file(filename,exp_type='ODNP_NMR_comp/field_dependent',
            expno=nodename,postproc=postproc,lookup=postproc_dict,fl=fl)
    s = s['t2':freq_slice]
    if s.get_prop('acq_params')['nPhaseSteps'] == 8:
        s.mean('nScans')
        s = s['ph1',1]['ph0',0].C
    else:
        s=s['ph1',1]['power',0]['nScans',0].C
    s.ift('t2')
    s.ft('t2')
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).C.convolve('t2',0.0001).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    s.ft('t2')
    fl.next('line plots')
    for z in range(len(s.getaxis('Field'))):
        fl.plot(abs(s['Field',z]),label='%d'%z)
    s_ = s['t2':field_slice].sum('t2')
    fl.next('sweep, without hermitian')
    fl.plot(abs(s_),'o-')
fl.show()
