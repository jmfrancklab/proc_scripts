from pyspecdata import *
from Utility import dBm2power
#to use type s = load_data("nameoffile")
def load_data(searchstr,exptype,expno):
    if (exptype=='ODNP'):
        files = search_filename(searchstr, 'test_equip')
        assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
        dirname, filename = os.path.split(files[0])
        nodename = 'signal'
        s = nddata_hdf5(filename+'/signal',
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
    if (exptype=='CPMG'):
        files = search_filename(searchstr, 'test_equip')
        assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
        dirname, filename = os.path.split(files[0])
        nodename = 'signal'
        s = nddata_hdf5(filename+'/signal',
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
    if (exptype=='Hahn_echoph'):
        files = search_filename(searchstr, 'test_equip')
        assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
        dirname, filename = os.path.split(files[0])
        nodename = 'signal'
        s = nddata_hdf5(filename+'/signal',
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
    if (exptype=='nutation'):
        files = search_filename(searchstr, 'test_equip')
        assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
        dirname, filename = os.path.split(files[0])
        nodename = 'nutation'
        s = nddata_hdf5(filename+'/nutation',
                directory=dirname)
        orig_t = s.getaxis('t')
        s.set_units('p_90','s')
        s.reorder('t',first=True)
        s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
        s.setaxis('ph2',r_[0.,2.]/4)
        s.setaxis('ph1',r_[0.,1.,2.,3.]/4)
        s.reorder('t2',first=False)
        return s
    if (exptype=='capture'):
        files = search_filename(searchstr, 'test_equip')
        assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
        dirname, filename = os.path.split(files[0])
        nodename = 'capture1'
        s = nddata_hdf5(filename+'/capture1',
                directory=dirname)
        s.set_units('t','s').name('Amplitude').set_units('V')
        return s
    if (exptype=='nmr'):
        files = search_filename(searchstr, 'test_equip')
        s = find_file(searchstr, exp_type='test_equip', dimname = 'indirect', expno=expno)
        return s
    if (exptype=='nmr_noecho'):
        files = search_filename(searchstr, 'test_equip')
        s = find_file(searchstr, exp_type='test_equip', dimname = 'indirect', expno=expno)
        fl = figlist_var()
        print(ndshape(s))
        print(s.getaxis('indirect'))
        s.chunk('indirect',['indirect','ph1','ph2'],[-1,4,2]) #expands the indirect dimension into indirect, ph1, and ph2. inner most dimension is the inner most in the loop in pulse sequence, is the one on the farthest right. brackets with numbers are the number of phase cycle steps in each one. the number of steps is unknown in 'indirect' and is therefore -1.
        print(s.getaxis('indirect'))
        print(s.getaxis('ph1'))
        print(s.getaxis('ph2'))
        s.setaxis('ph1',r_[0:4.]/4) #setting values of axis ph1 to line up
        s.setaxis('ph2',r_[0:2.]/4) #setting values of axis ph1 to line up
        s.setaxis('indirect', s.get_prop('vd'))
        fl.next('phased coherence domain') #switch to time domain as a string based name for fig
        s.ft(['ph1','ph2']) #fourier transforming from phase cycle dim to coherence dimension
        s.reorder(['indirect','t2'], first=False)
        fl.image(s) #labeling
        return s
