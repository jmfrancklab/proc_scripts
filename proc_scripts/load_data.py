"""Used to preprocess data based on type of experiment performed. Returns all data FTed
into the frequency domain with the phase cycles also FTed (coherence domain). Data is
not sliced or altered in anyway.

Parameters
==========
self:   nddata or h5 file

Returns
=======
nddata that has been FTed and in coherence domain
"""
from pyspecdata import *
from .Utility import dBm2power
import os
from sympy import symbols
import logging
import numpy as np
import logging
from pylab import *
from .DCCT_func import DCCT

#to use type s = load_data("nameoffile")
def proc_bruker_deut_IR_withecho_mancyc(s,fl=None):
    logging.info(strm("this is the 90 time"))
    if fl is not None:
        fl.next('raw data')
        fl.image(s.C.setaxis(
            'indirect','#').set_units('indirect','scan #'))
    s.chunk('indirect',['ph2','ph1','indirect'],[4,2,-1]) #expands the indirect dimension into indirect, ph1, and ph2. inner most dimension is the inner most in the loop in pulse sequence, is the one on the farthest right. Brackets with numbers are the number of phase cycle steps in each one. the number of steps is unknown in 'indirect' and is therefore -1.
    s.setaxis('ph1',r_[0:2.]/4) #setting values of axis ph1 to line up
    s.setaxis('ph2',r_[0:4.]/4) #setting values of axis ph1 to line up
    s.setaxis('indirect', s.get_prop('vd'))
    s.ft('t2',shift=True) #fourier transform
    if fl is not None:
        fl.next('IR prior to FTing ph')
        fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'))
    s.ft(['ph1','ph2']) #fourier transforming from phase cycle dim to coherence dimension
    s.reorder(['indirect','t2'], first=False)
    if fl is not None:
        s_forplot = s.C
        fl.next('FT')
        fl.image(s_forplot.C.setaxis('indirect','#').set_units('indirect','scan #'))
        fl.next('time domain (all $\\Delta p$)')
        s_forplot.ift('t2')
        fl.image(s_forplot.C.setaxis('indirect','#').set_units('indirect','scan #'))
        fl.next('frequency domain (all $\\Delta p$)')
        s_forplot.ft('t2',pad=4096)
        fl.image(s_forplot.C.setaxis('indirect','#').set_units('indirect','scan #'))
    return s

def proc_bruker_deut_IR_mancyc(s, fl=None):
    logging.info(strm("this is the d1",s.get_prop('acq')['D'][1]))
    if fl is not None:
        fl.next('raw data')
        fl.image(s)
    s.chunk('indirect',['indirect','ph1','ph2'],[-1,2,4]) #expands the indirect dimension into indirect, ph1, and ph2. 
                                                          #inner most dimension is the inner most in the loop in pulse sequence, 
                                                          #is the one on the farthest right. Brackets with numbers are the number 
                                                          #of phase cycle steps in each one. the number of steps is unknown in 
                                                          #'indirect' and is therefore -1.
    s.setaxis('ph1',r_[0:2.]/4) #setting values of axis ph1 to line up
    s.setaxis('ph2',r_[0:4.]/4) #setting values of axis ph1 to line up
    s.setaxis('indirect', s.get_prop('vd'))
#titling to coherence domain
    s.ft('t2',shift=True) #fourier transform
    s.ft(['ph1','ph2']) #fourier transforming from phase cycle dim to coherence dimension
    s.reorder(['indirect','t2'], first=False)
    if fl is not None:
        s_forplot = s.C
        s_forplot.setaxis('indirect','#').set_units('indirect','scan #')
        fl.next('FT + coherence domain')
        fl.image(s_forplot)
    if fl is not None:    
        fl.next('time domain (all $\\Delta p$)')
        s_forplot.ift('t2')
        fl.image(s_forplot)
    if fl is not None:    
        fl.next('frequency domain (all $\\Delta p$)')
        s_forplot.ft('t2',pad=4096)
        fl.image(s_forplot)
    return s

def proc_spincore_CPMG_v1(s, fl=None):
    logging.info("loading pre-processing for CPMG preprocessing")
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    nScans = s.get_prop('acq_params')['nScans']
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    deadtime_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank_s = s.get_prop('acq_params')['deblank_us']*1e-6
    marker_s = s.get_prop('acq_params')['marker_us']*1e-6
    tau1_s = s.get_prop('acq_params')['tau1_us']*1e-6
    pad_start_s = s.get_prop('acq_params')['pad_start_us']*1e-6
    pad_end_s = s.get_prop('acq_params')['pad_end_us']*1e-6
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    s.set_units('t','s')
    twice_tau = deblank_s + 2*p90_s + deadtime_s + pad_start_s + acq_time_s + pad_end_s + marker_s
    t2_axis = np.linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*twice_tau
    s.setaxis('nScans',r_[0:nScans])
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    # OK, so I made some changes and then realized that we need to assume that
    # the pulse sequence correctly balances the evolution between 2*p90_s/pi
    # (cavanagh chpt 3 this is the evolution during the 90 -- I'm not positive
    # if my expression is correct or not -- please do check/change, and leave
    # this comment in some form) and the center of the 180 pulse appropriately
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis)
    s.reorder(['ph1','tE','nScans','t2'])
    s.ft(['ph1'])
    s.reorder('nScans',first=True)
    if fl is not None:
        fl.next('time domain')
        fl.image(s)
    s.ft('t2', shift=True)
    if fl is not None:
        fl.next('frequency domain coh domain')
        fl.image(s)
    return s

def proc_bruker_T1CPMG_v1(s, fl=None):
    assert s.get_prop('acq')['L'][21] == 2, "phase cycle isn't correct!"
    assert s.get_prop('acq')['L'][22] == 4, "phase cycle isn't correct!"
    s.chunk('indirect',['indirect','ph1','ph2'],[-1,2,4])
    s.setaxis('ph1',r_[0,2]/4).setaxis('ph2',r_[0:4]/4)
    s.ft(['ph1','ph2'])
    s.reorder(['ph1','ph2','indirect'])
    if fl is not None:
        fl.next('raw data(t2,coh)')
        fl.image(s)
    #{{{removes CP aspect
    s.ift(['ph1','ph2'])
    s = s['ph2',[1,3]]
    #}}}
    s.setaxis('indirect', s.get_prop('vd'))
    s.reorder(['ph1','ph2','indirect','t2'])
    if fl is not None:
        fl.next('raw data with indirect set')
        fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'))
    anavpt_info = [j for j in s.get_prop('pulprog').split('\n') if 'anavpt' in j.lower()]
    anavpt_re = re.compile(r'.*\banavpt *= *([0-9]+)')
    anavpt_matches = (anavpt_re.match(j) for j in anavpt_info)
    for m in anavpt_matches:
        if m is not None:
            anavpt = int(m.groups()[0])
    actual_SW = 20e6/anavpt # JF: check that this is based on the manual's definition of anavpt
    bruker_final_t2_value = np.double(s.getaxis('t2')[-1].item())
    s.setaxis('t2',1./actual_SW*r_[0:ndshape(s)['t2']]) # reset t2 axis to true values based on anavpt
    logging.debug(strm("the final t2 value according to the Bruker SW_h was",
        bruker_final_t2_value, "but I determine it to be",
        np.double(s.getaxis('t2')[-1].item()), "with anavpt"))
    nEchoes = s.get_prop('acq')['L'][25]
    dwdel1 = s.get_prop('acq')['DE']*1e-6
    dwdel2 = (anavpt*0.05e-6)/2
    #d12 is read as 0 if taken from parameters bc its too small
    d12 = s.get_prop('acq')['D'][12]
    d11 = s.get_prop('acq')['D'][11]
    p90_s = s.get_prop('acq')['P'][1]*1e-6
    quad_pts = ndshape(s)['t2'] # note tha twe have not yet chunked t2
    nPoints = quad_pts/nEchoes
    acq_time = dwdel2*nPoints*2
    # {{{ these are hard-coded for the pulse sequence
    #     if we need to, we could pull these from the pulse sequence, as we do
    #     for anavpt above
    tau_extra = d12
    tau_pad_start = tau_extra-dwdel1-6e-6
    tau_pad_end = tau_extra-6e-6
    twice_tau = 2*p90_s + 5e-6 + tau_pad_start + 1e-6 + acq_time + tau_pad_end +1e-6
    # twice_tau should be the period from one 180 to another
    # }}}
    s.chunk('t2',['tE','t2'],[nEchoes,-1])
    s.setaxis('tE', (1+r_[0:nEchoes])*twice_tau)
    s.ft('t2', shift=True)
    s.ft(['ph1','ph2'])
    s.reorder(['ph1','ph2','indirect'])
    if fl is not None:
        fl.next('freq domain coh domain')
        fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'))
    s.ift('t2')
    if fl is not None:
        fl.next('t2 chunked', figsize=(5,20))
        fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'))
    s.ft('t2')
    return s

def proc_bruker_CPMG_v1(s,fl=None):
    s.chunk('indirect',['ph1','ph2','indirect'],[4,2,-1])
    s.setaxis('ph1',r_[0:4]/4.)
    s.setaxis('ph2',r_[0:2]/2.)
    if fl is not None:
        fl.next('raw data before')
        fl.image(s)
    s.ft(['ph1','ph2'])
    s.reorder(['ph1','ph2','indirect','t2'])
    anavpt_info = [j for j in s.get_prop('pulprog').split('\n') if 'anavpt' in j.lower()]
    anavpt_re = re.compile(r'.*\banavpt *= *([0-9]+)')
    anavpt_matches = (anavpt_re.match(j) for j in anavpt_info)
    for m in anavpt_matches:
        if m is not None:
            anavpt = int(m.groups()[0])
    actual_SW = 20e6/anavpt # JF: check that this is based on the manual's definition of anavpt
    bruker_final_t2_value = np.double(s.getaxis('t2')[-1].item())
    s.setaxis('t2',1./actual_SW*r_[0:ndshape(s)['t2']]) # reset t2 axis to true values based on anavpt
    logging.debug(strm("the final t2 value according to the Bruker SW_h was",
        bruker_final_t2_value, "but I determine it to be",
        np.double(s.getaxis('t2')[-1].item()), "with anavpt"))
    nEchoes = s.get_prop('acq')['L'][25]
    dwdel1 = s.get_prop('acq')['DE']*1e-6
    dwdel2 = (anavpt*0.05e-6)/2
    #d12 is read as 0 if taken from parameters bc its too small
    d12 = 20e-6     
    d11 = s.get_prop('acq')['D'][11]
    p90_s = s.get_prop('acq')['P'][1]*1e-6
    quad_pts = ndshape(s)['t2'] # note tha twe have not yet chunked t2
    nPoints = quad_pts/nEchoes
    acq_time = dwdel2*nPoints*2
    # {{{ these are hard-coded for the pulse sequence
    #     if we need to, we could pull these from the pulse sequence, as we do
    #     for anavpt above
    tau_extra = 20e-6
    tau_pad_start = tau_extra-dwdel1-6e-6
    tau_pad_end = tau_extra-6e-6
    twice_tau = 2*p90_s + 5e-6 + tau_pad_start + 1e-6 + acq_time + tau_pad_end +1e-6
    # twice_tau should be the period from one 180 to another
    # }}}
    s.set_units('t2','us')
    s.chunk('t2',['tE','t2'],[nEchoes,-1])
    s.setaxis('tE', (1+r_[0:nEchoes])*twice_tau)
    s.ft('t2', shift=True)
    s.reorder(['ph1','ph2','indirect'])
    if fl is not None:
        fl.next('freq domain coh domain')
        fl.image(s)
    s.ift('t2')
    if fl is not None:
        fl.next('t2 chunked', figsize=(5,20))
        fl.image(s)
    s.ft('t2')
    return s

def proc_Hahn_echoph(s, fl=None):
    logging.info("loading pre-processing for Hahn_echoph")
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = 8 
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph2','ph1'])
    s.setaxis('nScans',r_[0:nScans])
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    if fl is not None:
        fl.next('raw data, chunked')
        fl.image(abs(s))
    s.ft(['ph1','ph2'])
    if fl is not None:
        fl.next('coherence')
        fl.image(abs(s))
    return s

def proc_spincore_IR(s,fl=None):
    vd_axis = s.getaxis('vd')
    s.reorder(['ph2','ph1']).set_units('t2','s')
    s.ft('t2', shift=True)
    s.ft(['ph2','ph1'])
    if fl is not None:
        fl.next('raw data -- coherence channels')
        fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    s.ift('t2')
    if fl is not None:
        fl.next('time domain (all $\\Delta p$)')
        fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'))
    s.ft('t2', pad=4096)
    if fl is not None:
        fl.next('frequency domain (all $\\Delta p$)')
        fl.image(s.C.setaxis('vd','#').set_units('vd','scan #'),black=False)
    return s

def proc_nutation(s,fl=None):
    logging.info("loading pre-processing for nutation")
    orig_t = s.getaxis('t')
    s.set_units('p_90','s')
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,2,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,2.]/4)
    s.reorder('t2',first=False)
    s.ft(['ph2','ph1'])
    if fl is not None:
        fl.next('after phase cycle FT')
        fl.image(s['ph1',1]['ph2',0].C.human_units())
    s.ft('t2',shift=True)
    if fl is not None:
        fl.next('freq domain')
        fl.image(s)
    return s

def proc_nutation_amp(s,fl=None):
    logging.info("loading pre-processing for nutation")
    s.set_units('t2','s')
    s.ft('t2',shift=True)
    if fl is not None:
        fl.next('look for drift')
        fl.image(s.C.smoosh(['ph2','ph1'],'transient').reorder('transient').setaxis('transient','#').run(abs),
                interpolation='bilinear')
    s.reorder(['ph1','ph2'])
    s.setaxis('ph2',r_[0:2]/4).setaxis('ph1',r_[0:4]/4)
    if 'p_90' in s.dimlabels:
        s.set_units('p_90','s')
    s.ft(['ph1','ph2'])
    return s

def proc_nutation_chunked(s,fl=None):
    logging.info("loading pre-processing for nutation")
    s.reorder(['ph1','ph2'])
    s.set_units('t2','s')
    s.set_units('p_90','s')
    #s.reorder('t2',first=True)
    s.ft(['ph1','ph2'])
    s.reorder(['ph1','ph2','p_90'])
    if fl is not None:
        fl.next('Raw Data - Time Domain')
        fl.image(s.C.human_units())
    s.ft('t2',shift=False)
    if fl is not None:
        fl.next('Raw Data- Frequency Domain')
        fl.image(s)
    return s

def proc_var_tau(s,fl=None):
    s.get_prop('SW')
    if 'ph1' not in s.dimlabels:
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
        s.setaxis('ph2',r_[0,2]/4)
        s.setaxis('ph1',r_[0:4]/4)
    s.set_units('t2','s') # this should already be set -- why not?
    s *= 2e-6/1.11e4 # convert from SpinCore to V (amp)
    s.set_units('V')
    if fl is not None:
        fl.next('raw signal!')
    s.ft('t2', shift=True).ft(['ph1','ph2'])
    s.reorder(['ph1','ph2','tau'])
    if fl is not None:
        fl.plot(abs(s).smoosh(['ph2','ph1','tau'],'transients'), alpha=0.2)
        fl.next('raw signal')
        fl.image(s)
    return s

def proc_spincore_ODNP_v1(s,fl=None):
    logging.info("loading pre-processing for ODNP")
    prog_power = s.getaxis('power').copy()
    logging.info(strm("programmed powers",prog_power))
    s.setaxis('power',r_[
        0,dBm2power(np.array(s.get_prop('meter_powers'))+20)]
        ).set_units('power','W')
    logging.info(strm("meter powers",s.get_prop('meter_powers')))
    logging.info(strm("actual powers",s.getaxis('power')))
    logging.info(strm("ratio of actual to programmed power",
               s.getaxis('power')/prog_power))
    nPoints = s.get_prop('acq_params')['nPoints']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    s.set_units('t','s')
    s.chunk('t',['ph2','ph1','t2'],[2,2,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,2.]/4})
    s.reorder(['ph2','ph1'])
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2']) # Fourier Transforms coherence channels
    if fl is not None:
        fl.next('all data: frequency domain')
        fl.image(s.C.setaxis('power','#').set_units('power','scan #'))
    return s

def proc_capture(s):
    logging.info("loading pre-processing for square wave capture")
    s.set_units('t','s').name('Amplitude').set_units('V')
    return s

def proc_DOSY_CPMG(s):
    logging.info("loading pre-processing for DOSY-CPMG")
    # {{{ all of this would be your "preprocessing" and would be tied to the name of your pulse sequence
    l22 = int(s.get_prop('acq')['L'][22]) # b/c the l are integers by definition
    l25 = int(s.get_prop('acq')['L'][25])
    d12 = s.get_prop('acq')['D'][12]
    d11 = s.get_prop('acq')['D'][11]
    p1 = s.get_prop('acq')['P'][1]
    ppg = s.get_prop('pulprog')
    # {{{ these are explanatory -- maybe comment them out?
    m = re.search(('.*dwdel1=.*'),ppg,flags=re.IGNORECASE)
    logging.info(strm(m.groups())) # show the line that sets dwdel1
    # then look for de and depa
    logging.info(strm([(j,s.get_prop('acq')[j]) for j in s.get_prop('acq').keys() if 'de' in j.lower()]))
    # I actually can't find depa
    # }}}
    m = re.search('\ndefine list<grad_scalar> gl1 = {(.*)}',ppg)
    grad_list = array([float(j.group()) for j in re.finditer('([0-9.]+)',m.groups()[0])])
    m = re.search('([0-9.]+) G/mm', s.get_prop('gradient_calib'))
    grad_list *= float(m.groups()[0])*0.1
    dwdel1 = 3.5e-6 # where does this come from? DE is actually larger than this?
    # {{{ find anavpt without hard-setting
    m = re.search('"anavpt=([0-9]+)"',ppg)
    if m is None:
        raise ValueError("I can't find anavpt in the pulse sequence")
    anavpt = int(m.groups()[0])
    # }}}
    dwdel2 = (anavpt*0.05e-6)/2
    TD = s.get_prop('acq')['TD2']
    quadrature_points = TD/2
    num_points_per_echo = quadrature_points/l25
    acq_time = dwdel2*num_points_per_echo*2
    # {{{ so, in principle, later, we can/should do what I did above (w/ eval),
    # but it's getting crazy now, so I stop for now
    tau_extra = 20e-6
    tau_pad = tau_extra-6e-6
    tau_pad_start = tau_extra-dwdel1-6e-6
    tau_pad_end = tau_extra-6e-6
    tE = dwdel1 + 5e-6 + tau_pad_start + 1e-6 + num_points_per_echo*(dwdel2*2) + tau_pad_end
    # }}}
    s.chunk('indirect',['indirect','phcyc'],[l22,-1])
    s.chunk('phcyc',['ph8','ph4','m','n'],[2,2,2,2])
    s.setaxis('ph8',r_[0.,2.]/4)
    s.setaxis('ph4',r_[0.,2.]/4)
    s.setaxis('m',r_[0,2.]/4)
    s.setaxis('n',r_[0,2.]/4)
    s.ft(['ph8','ph4','m','n'])
    s.reorder(['m','n','ph4','ph8','indirect','t2'])
    s.setaxis('indirect',grad_list)
    fl.next('abs raw data')
    fl.image(abs(s))
    s.chunk('t2',['echo','t2'],[l25,-1])
    s.reorder(['m','n','ph4','ph8','indirect','echo','t2'])
    s.ft('t2', shift=True).ift('t2') # this is overkill -- need a pyspecdata function that does this w/out the fft
    # }}}
    return s

def proc_ESR(s):
    logging.info("loading preprocessing for ESR linewidth calculation")
    s -= s['$B_0$',:50].C.mean('$B_0$')
    s_integral = s.C.run_nopop(np.cumsum,'$B_0$')
    x1,x2 = s_integral.getaxis('$B_0$')[r_[5,-5]]
    y1 = s_integral.data[:5].mean()
    y2 = s_integral.data[-5:].mean()
    straight_baseline = (s.fromaxis('$B_0$')-x1)*(y2-y1)/(x2-x1)
    s_integral -= straight_baseline
    s_integral /= s_integral.data.mean()
    center_field = (s_integral * s.fromaxis('$B_0$')).mean('$B_0$').item()
    s.setaxis('$B_0$',lambda x: x-center_field)
    s_integral = s.C.run_nopop(np.cumsum,'$B_0$')
    logging.info(strm(s_integral))
    return s    

postproc_dict = {'ag_IR2H':proc_bruker_deut_IR_withecho_mancyc,
        'ab_ir2h':proc_bruker_deut_IR_mancyc,
        'ag_CPMG_strob':proc_bruker_CPMG_v1,
        'ag_T1CPMG_2h':proc_bruker_T1CPMG_v1,
        'chirp':proc_capture,
        'spincore_CPMG_v1':proc_spincore_CPMG_v1,
        'spincore_Hahn_echoph_v1':proc_Hahn_echoph,
        'spincore_IR_v1':proc_spincore_IR,
        'spincore_nutation_v1':proc_nutation,
        'spincore_nutation_v2':proc_nutation_amp,
        'spincore_nutation_v3':proc_nutation_chunked,
        'spincore_ODNP_v1':proc_spincore_ODNP_v1,
        'spincore_var_tau_v1':proc_var_tau,
        'square_wave_capture_v1':proc_capture,
        'DOSY_CPMG_v1':proc_DOSY_CPMG,
        'ESR_linewidth':proc_ESR,

}
