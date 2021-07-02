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
filename = '210702_500uM_TEMPO_hexane_cap_probe_field_dep' #'210611_S175R1a_pR_DDM_field_dep'
gamma_eff = (14.824903/3489.4)#(14.893851/3505.6) # MHz / G
f_dip = 9.8214286#9.82103 # GHz
for nodename,postproc,label_str,freq_slice,field_slice in [
        ('field_sweep_2',#'32dBm_finer',
        'field_sweep','TEMPO field sweep',(-500,700),(-200,500)),
        ]:
    s = find_file(filename,exp_type='ODNP_NMR_comp/field_dependent',
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
    #fl.show();quit()
    s_ = s['t2':field_slice].sum('t2')
    fl.next('sweep, without hermitian')
    fl.plot(abs(s_),'o-')
    field_idx = (abs(s_.data)).argmax()
    print('At $B_0$ = %0.1f, $f_0$ = %0.8e'%(s.getaxis('Field')[field_idx],freqs[field_idx]))
    fitting = abs(s_).polyfit('Field',order=2)
    Field = nddata(r_[3503.5:3508:100j],'Field')
    #fl.plot(Field.eval_poly(fitting,'Field'),label='fit')
    print('I found a max at',Field.eval_poly(fitting,'Field').argmax().item())
fl.show()
