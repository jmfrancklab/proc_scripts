from pyspecdata import *
from Utility import dBm2power
#to use type s = load_data("nameoffile")
def load_data(searchstr):
    files = search_filename(searchstr, 'test_equip')
    assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
    dirname, filename = os.path.split(files[0])
    nodename = 'signal'
    s = nddata_hdf5(filename+'/signal',
            directory=dirname)
    #return s
    s.get_prop()
    if 'meter_powers' in s.get_prop():
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
    else:
        return s
