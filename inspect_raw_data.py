from pyspecdata import *
import os
from sympy import symbols
import matplotlib.pyplot as plt
import numpy as np
t2 = symbols('t2')
fig_dict={}
clock_correction = 1.0829/998.253
for j in ['raw data','select coherence']:
    fig,axlist = plt.subplots(nrows=3,ncols=2)
    fig_dict[j] = axlist.ravel()
for j in range(1,7):
    fname = search_filename("200122_IR_water_%s"%j,'test_equip')
    assert len(fname)==1
    fname = os.path.split(os.path.normpath(fname[0]))
    dirname,fname = fname
    logger.info(strm("looking for nutation node in",fname,"in",dirname))
    s = nddata_hdf5(fname+'/signal', directory=dirname)
    # data should be chunked before storage
    s.chunk('t',['ph2','ph1','t2'],[4,2,-1])
    s.setaxis('ph1',r_[0:4:2]/4.)
    s.setaxis('ph2',r_[0:4]/4.)
    s.set_units('t2','s')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    s.ft('t2',shift=True).reorder(['ph2','ph1'])
    s *= np.exp(-1j*s.fromaxis('vd')*clock_correction)
    plt.sca(fig_dict['raw data'][j-1])
    image(s['t2':(-200,200)].C.setaxis('vd','#').set_units('vd','scan #'),
            interpolation='bicubic')
    plt.sca(fig_dict['select coherence'][j-1])
    s = s['ph1',0]['ph2',1]
    image(s['t2':(-200,200)].C.setaxis('vd','#').set_units('vd','scan #'))
show()
