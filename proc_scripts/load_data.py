from pyspecdata import *
from .Utility import dBm2power
from proc_scripts import *
import os
from sympy import symbols
import logging
fl=figlist_var()
#to use type s = load_data("nameoffile")
def proc_bruker_90_pulse(s,fl=None):
    s.ft('t2',shift=True)
    if fl is not None:
        fl.next('raw-frequency domain')
        fl.plot(s)
    if fl is not None:
        fl.next('raw-time domain')
        s.ift('t2')
        fl.plot(s)
    return s
def proc_bruker_deut_IR_withecho_mancyc(s,fl=fl):
    logger.info(strm("this is the 90 time",s.get_prop('acq')['P'][1]))
    s.chunk('indirect',['ph2','ph1','indirect'],[2,4,-1]) #expands the indirect dimension into indirect, ph1, and ph2. inner most dimension is the inner most in the loop in pulse sequence, is the one on the farthest right. Brackets with numbers are the number of phase cycle steps in each one. the number of steps is unknown in 'indirect' and is therefore -1.
    s.setaxis('ph1',r_[0:4.]/4) #setting values of axis ph1 to line up
    s.setaxis('ph2',r_[0:2.]/4) #setting values of axis ph1 to line up
    s.setaxis('indirect', s.get_prop('vd'))
#titling to coherence domain
    s.ft('t2',shift=True) #fourier transform
    if fl is not None:
        fl.next('IR prior to FTing ph')
        fl.image(s)
    s.ft(['ph1','ph2'])
    s.reorder(['indirect','t2'], first=False)
    if fl is not None:
        s_forplot = s.C
        fl.next('FT')
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

def proc_bruker_deut_IR_mancyc(s, fl=None):
    s.chunk('indirect',['indirect','ph1','ph2'],[-1,4,2]) #expands the indirect dimension into indirect, ph1, and ph2. inner most dimension is the inner most in the loop in pulse sequence, is the one on the farthest right. Brackets with numbers are the number of phase cycle steps in each one. the number of steps is unknown in 'indirect' and is therefore -1.
    s.setaxis('ph1',r_[0:4.]/4) #setting values of axis ph1 to line up
    s.setaxis('ph2',r_[0:2.]/4) #setting values of axis ph1 to line up
    s.setaxis('indirect', s.get_prop('vd'))
#titling to coherence domain
    s.ft('t2',shift=True) #fourier transform
    s.ft(['ph1','ph2']) #fourier transforming from phase cycle dim to coherence dimension
    s.reorder(['indirect','t2'], first=False)
    if fl is not None:
        s_forplot = s.C
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
    logger.info("loading pre-processing for CPMG preprocessing")
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
    t2_axis = linspace(0,acq_time_s,nPoints)
    s.setaxis('nScans',r_[0:nScans])
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    # OK, so I made some changes and then realized that we need to assume that
    # the pulse sequence correctly balances the evolution between 2*p90_s/pi
    # (cavanagh chpt 3 this is the evolution during the 90 -- I'm not positive
    # if my expression is correct or not -- please do check/change, and leave
    # this comment in some form) and the center of the 180 pulse appropriately
    s.setaxis('tE', (1+r_[0:nEchoes])*twice_tau)
    s.setaxis('ph1',r_[0.,2.]/4)
    s.ft('t2', shift=True)
    if fl is not None:
        fl.next('raw data - chunking ft')
        fl.image(s)
    s.ft(['ph1'])
    s.ift('t2')
    s.reorder('nScans',first=True)
    s = s['ph1',1].C
    s.mean('nScans')
    s.reorder('t2',first=True)
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

def proc_spincore_IR(s,clock_correction,fl=None):
    s.rename('vd','indirect')
    s.reorder(['ph1','ph2','indirect','t2'])
    s.ft(['ph2','ph1'])
    s.ft('t2', shift=True)
    if fl is not None:
        fl.next('raw data -- coherence channels')
        fl.image(s)
    s.ift('t2')
    if fl is not None:
        fl.next('time domain (all $\\Delta p$)')
        fl.image(s)
    s.ft('t2', pad=4096)
    if fl is not None:
        fl.next('frequency domain (all $\\Delta p$)')
        fl.image(s)
    s *= exp(-1j*s.fromaxis('indirect')*clock_correction)
    return s

def proc_nutation(s):
    logging.info("loading pre-processing for nutation")
    orig_t = s.getaxis('t')
    s.set_units('p_90','s')
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.reorder('t2',first=False)
    s.ft('t2',shift=True)
    s.ift('t2')
    return s

def proc_spincore_ODNP_v1(s,fl=None):
    logging.info("loading pre-processing for ODNP")
    prog_power = s.getaxis('power').copy()
    logging.info(strm("programmed powers",prog_power))
    s.setaxis('power',r_[
        0,dBm2power(array(s.get_prop('meter_powers'))+20)]
        ).set_units('power','W')
    logging.info(strm("meter powers",s.get_prop('meter_powers')))
    logging.info(strm("actual powers",s.getaxis('power')))
    logging.info(strm("ratio of actual to programmed power",
               s.getaxis('power')/prog_power))
    nPoints = s.get_prop('acq_params')['nPoints']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    s.set_units('t','s')
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph2','ph1'])
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2']) # Fourier Transforms coherence channels
    if fl is not None:
        fl.next('all data: frequency domain')
        fl.image(s)
    return s

def proc_square_wave_capture(s):
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

postproc_dict = {'zg2h':proc_bruker_90_pulse,
        'ag_IR2H':proc_bruker_deut_IR_withecho_mancyc,
        'ab_ir2h':proc_bruker_deut_IR_mancyc,
        'spincore_CPMG_v1':proc_spincore_CPMG_v1,
        'spincore_Hahn_echoph_v1':proc_Hahn_echoph,
        'spincore_nutation_v1':proc_nutation,
        'spincore_IR_v1':proc_spincore_IR,
        'spincore_ODNP_v1':proc_spincore_ODNP_v1,
        'square_wave_capture_v1':proc_square_wave_capture,
        'DOSY_CPMG_v1':proc_DOSY_CPMG}

