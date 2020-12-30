from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from pylab import *
import h5py as h5
fl = figlist_var()
GDS = True
for filename,expno,f_range in [
        ('201229_input_referred_3','signal',(14.7e3,15e3))
        ]:
    fl.basename=filename
    print('analyizing',filename)
    fullname = search_filename(filename,
            exp_type='test_equip',
            unique=True)
    with h5.File(fullname,'r') as fp:
        print(fp.keys())
    s = find_file(filename,
            exp_type='test_equip',expno=expno)
    if GDS:
        s.set_units('t','s')
        s.set_units('V')
        print(ndshape(s))
        fl.next('show signal!')
        s.ft('t',shift=True)
        s = s['t':(-34e3,34e3)]
        ax1 = subplot(2,1,1)
        fl.plot(abs(s),alpha=0.2)
        #fl.show();quit()
        ax2=subplot(2,1,2)
        fl.plot(abs(s)['t':(14e3,15e3)],alpha=0.2)
        s['t':(None,f_range[0])] = 0
        s['t':(f_range[1],None)] = 0
        s *= 2
        fl.plot(abs(s), ':',ax=ax1)
        fl.plot(abs(s)['t':(14e3,15e3)],alpha=0.2,ax=ax2)
        fl.next('time domain')
        s.setaxis('t', lambda x: x-14e6)
        s.ift('t')
        fl.plot(abs(s),alpha=0.5, linewidth=2)
        fl.plot(s.real, alpha = 0.2, linewidth = 0.5)
        fl.plot(s.imag, ':', alpha=0.2, linewidth = 0.5)
    else:
        print(ndshape(s))
        s.get_prop('SW')
        if 'ph1' not in s.dimlabels:
            s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
            s.setaxis('ph2',r_[0,2]/4)
            s.setaxis('ph1',r_[0:4]/4)
        print(ndshape(s))
        s.set_units('t2','s')
        s.set_units('Spincore')
        fl.next('raw signal!')
        s.ft('t2',shift=True).ft(['ph1','ph2'])
        s.reorder(['ph1','ph2'])#,'tau'])
        fl.plot(abs(s).smoosh(['ph2','ph1'],'transients'),alpha=0.2)
        fl.next('raw signal')
        fl.image(s)
        #fl.show();quit()
        s = s['t2':f_range]
        s = s['ph1',+1]['ph2',-2]
        s.ift('t2')
        fl.next('echoes')
        fl.plot(s.real,alpha=0.2,linewidth=0.5)
        fl.plot(abs(s),alpha=0.5,linewidth=1)
        fl.plot(1.3e3, 'x',label='estimated signal at 0')
    fl.show();quit()
