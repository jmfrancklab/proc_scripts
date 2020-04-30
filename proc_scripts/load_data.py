from pyspecdata import *
from Utility import dBm2power
#to use type s = load_data("nameoffile")
def proc_bruker_deut_IR_withecho_mancyc(s):
    raise RuntimeError("this is where postprocessing would be implemented -- not implemented yet")
def proc_bruker_deut_IR_mancyc(s):
    raise RuntimeError("this is where postprocessing would be implemented -- not implemented yet")
postproc_lookup = {'ag_IR2H':proc_bruker_deut_IR_withecho_mancyc,
        'ab_ir2h':proc_bruker_deut_IR_mancyc}
def load_data(searchstr, exp_type, which_exp=None, postproc=None):
    if postproc=='spincore_ODNP_v1':
        filename = search_filename(searchstr, exp_type, unique=True)
        dirname, filename = os.path.split(filename)
        nodename = 'signal'
        s = nddata_hdf5(filename+'/'+which_exp,
                directory=dirname)
        print("getting acquisition parameters")
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
        return s
    elif postproc=='CPMG':
        filename = search_filename(searchstr, exp_type)
        dirname, filename = os.path.split(filename)
        nodename = 'signal'
        s = nddata_hdf5(filename+'/'+which_exp,
                directory=dirname)
        print("getting acquisition parameters")
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
        tE_axis = r_[1:nEchoes+1]*twice_tau
        s.setaxis('t',None)
        s.setaxis('nScans',r_[0:nScans])
        s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
        s.setaxis('ph1',r_[0.,2.]/4)
        s.setaxis('tE',tE_axis)
        s.setaxis('t2',t2_axis)
        return s
    elif postproc=='Hahn_echoph':
        filename = search_filename(searchstr, exp_type)
        dirname, filename = os.path.split(filename)
        nodename = 'signal'
        s = nddata_hdf5(filename+'/'+which_exp,
                directory=dirname)
        print("getting acquisition parameters")
        nPoints = s.get_prop('acq_params')['nPoints']
        nEchoes = s.get_prop('acq_params')['nEchoes']
        nPhaseSteps = 8 
        SW_kHz = s.get_prop('acq_params')['SW_kHz']
        nScans = s.get_prop('acq_params')['nScans']
        print(ndshape(s))
        s.reorder('t',first=True)
        t2_axis = s.getaxis('t')[0:128]
        s.setaxis('t',None)
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
        s.setaxis('ph2',r_[0.,2.]/4)
        s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
        s.setaxis('t2',t2_axis)
        s.setaxis('nScans',r_[0:nScans])
        s.reorder('t2',first=False)
        return s
    elif postproc=='nutation':
        filename = search_filename(searchstr, exp_type)
        dirname, filename = os.path.split(filename)
        nodename = 'nutation'
        s = nddata_hdf5(filename+'/'+which_exp,
                directory=dirname)
        orig_t = s.getaxis('t')
        s.set_units('p_90','s')
        s.reorder('t',first=True)
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
        s.setaxis('ph2',r_[0.,2.]/4)
        s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
        s.reorder('t2',first=False)
        return s
    elif postproc=='capture':
        filename = search_filename(searchstr, exp_type)
        dirname, filename = os.path.split(filename)
        nodename = 'capture1'
        s = nddata_hdf5(filename+'/'+which_exp,
                directory=dirname)
        s.set_units('t','s').name('Amplitude').set_units('V')
        return s
    elif postproc is None:
        print("You left postproc unset, so I'm assuming you're going to let me choose what to do.  Right now, this only works for Bruker format files")
        # if we set s.set_prop('postproc_type'...), then find_file should automatically recognize what to do
        s = find_file(searchstr, exp_type=exp_type, dimname='indirect',
                expno=which_exp, postproc_lookup=postproc_lookup)
        return s
    else:
        raise ValueError("I can't determine the type of postprocessing")
