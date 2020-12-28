from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
for date,id_string in [
        ('201209','with_tune_limit_4uV_Vout'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'capture1'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip'))
    print(ndshape(s))
    fl.next('raw without tune limiter')
    fl.plot(s)
    s.ft('t')
    fl.next('FTed abs')
    fl.plot(abs(s))
    fl.show();quit()

