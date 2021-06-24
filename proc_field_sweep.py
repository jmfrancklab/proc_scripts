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
# sphinx_gallery_thumbnail_number = 2
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery

fl = figlist_var()
t2 = symbols('t2')
filter_bandwidth = 20e3
filename = '210615_S175R1a_pR_DDM_field_dep_2' #'210611_S175R1a_pR_DDM_field_dep'
gamma_eff = (14.921343/3512.1)#(14.893851/3505.6) # MHz / G
for nodename,postproc,label_str,freq_slice,field_slice in [
        ('32dBm_real_centered',#'32dBm_finer',
        'field_sweep','Sams field sweep',(-700,700),(-400,300)),
        ]:
    s = find_file(filename,exp_type='odnp',#'ODNP_NMR_comp/field_dependent',
            expno=nodename,postproc=postproc,lookup=postproc_dict,fl=fl)
    freqs = s.get_prop('acq_params')['mw_freqs']
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
    field_idx = (abs(s_.data)).argmax()
    print(field_idx)
    print('At %0.1f G, the NMR frequency is %0.8e MHz'%(s.getaxis('Field')[field_idx],freqs[field_idx]))
#    gamma_eff = (/s.getaxis('Field')[field_idx])
#    s_.setaxis('Field',lambda x: gamma_eff*x)
#    fl.next('sweep across carrier freq')
#    fl.plot(abs(s_),'go-')
#    xlabel('frequency (MHz)')
#    s_.setaxis('Field',lambda x: x/f_dip)
#    fl.next('sweep across ppt-value')
#    fl.plot(abs(s_),'mo-')
#    xlabel('ppt (MHz/GHz)')
fl.show()
