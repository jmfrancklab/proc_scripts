from pyspecdata import *
from .Utility import dBm2power
fl = figlist_var()
#to use type s = load_data("nameoffile")
def proc_bruker_deut_IR_withecho_mancyc(s):
    raise RuntimeError("this is where postprocessing would be implemented -- not implemented yet")

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
        fl.next('time domain (all $\\Delta p$)')
        s_forplot.ift('t2')
        fl.image(s_forplot)
        fl.next('frequency domain (all $\\Delta p$)')
        s_forplot.ft('t2',pad=4096)
        fl.image(s_forplot)
    return s
    #raise RuntimeError("this is where postprocessing would be implemented -- not implemented yet")

def proc_spincore_CPMG_v1(s, fl=None):
    logging.basicConfig()
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
    print("loading pre-processing for Hahn_echoph")
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = 8 
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nScans = s.get_prop('acq_params')['nScans']
    s.reorder('t',first=True)
    t2_axis = s.getaxis('t')[0:256]
    s.setaxis('t',None)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph2','ph1'])
    s.setaxis('t2',t2_axis)
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

def proc_spincore_IR(s):
    s.rename('vd','indirect')
    s.reorder(['ph1','ph2','indirect','t2'])
    fl.next('raw data -- coherence channels')
    s.ft(['ph2','ph1'])
    s.ft('t2', shift=True)
    fl.image(s)
    fl.next('time domain (all $\\Delta p$)')
    s.ift('t2')
    fl.image(s)
    fl.next('frequency domain (all $\\Delta p$)')
    s.ft('t2', pad=4096)
    fl.image(s)
    return s

def proc_nutation(s):
    print("loading pre-processing for nutation")
    orig_t = s.getaxis('t')
    s.set_units('p_90','s')
    s.reorder('t',first=True)
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.setaxis('ph2',r_[0.,2.]/4)
    s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
    s.reorder('t2',first=False)
    return s

def proc_spincore_ODNP_v1(s):
    print("loading pre-processing for ODNP")
    prog_power = s.getaxis('power').copy()
    print("programmed powers",prog_power)
    s.setaxis('power',r_[
        0,dBm2power(array(s.get_prop('meter_powers'))+20)]
        ).set_units('power','W')
    print("meter powers",s.get_prop('meter_powers'))
    print("actual powers",s.getaxis('power'))
    print("ratio of actual to programmed power",
               s.getaxis('power')/prog_power)
    nPoints = s.get_prop('acq_params')['nPoints']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    s.set_units('t','s')
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph2','ph1'])
    s.ft('t2',shift=True)
    return s

def proc_square_wave_capture(s):
    print("loading pre-processing for square wave capture")
    s.set_units('t','s').name('Amplitude').set_units('V')
    return s

def proc_90_pulse(s):
    return s

postproc_dict = {'ag_IR2H':proc_bruker_deut_IR_withecho_mancyc,
        'ab_ir2h':proc_bruker_deut_IR_mancyc,
        'spincore_CPMG_v1':proc_spincore_CPMG_v1,
        'spincore_Hahn_echoph_v1':proc_Hahn_echoph,
        'spincore_nutation_v1':proc_nutation,
        'spincore_IR_v1':proc_spincore_IR,
        'spincore_ODNP_v1':proc_spincore_ODNP_v1,
        'square_wave_capture_v1':proc_square_wave_capture,
        'zg2h':proc_90_pulse}

