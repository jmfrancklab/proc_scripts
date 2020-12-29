from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
fl = figlist_var()
for date,id_string in [
        ('201229','input_referred_3'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(exp_type = 'test_equip'))
    print(ndshape(s))
    fl.next('Spincore with tuned limiter')
    fl.plot(s)
    s.ft('t')
    #fl.next('FTed abs')
    #fl.plot(s)
    fl.show();quit()

