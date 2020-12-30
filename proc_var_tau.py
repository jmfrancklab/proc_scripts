from pylab import *
from pyspecdata import *
import h5py as h5

with figlist_var() as fl:
    for filename,expno,f_range in [
            ("201209_Ni_sol_probe_var_tau",'var_tau',(-13.5e3,0))
            ]:
        fl.basename = filename
        print("analyzing",filename)
        # {{{ I wanted to see what your dataset is called
        fullname = search_filename(filename,
                exp_type='var_tau',
                unique=True)
        with h5.File(fullname,'r') as fp:
            print(fp.keys())
        # }}}
        d = find_file(filename,
                exp_type='var_tau',
                expno=expno)
        print(ndshape(d))
        d.get_prop('SW')
        if 'ph1' not in d.dimlabels:
            d.chunk('t',['ph2','ph1','t2'],[2,4,-1])
            d.setaxis('ph2',r_[0,2]/4)
            d.setaxis('ph1',r_[0:4]/4)
        print(ndshape(d))
        d.set_units('t2','s') # this should already be set -- why not?
        d *= 2e-6/1.11e4 # convert from SpinCore to V (amp)
        d.set_units('V')
        fl.next('raw signal!')
        d.ft('t2', shift=True).ft(['ph1','ph2'])
        d.reorder(['ph1','ph2','tau'])
        fl.plot(abs(d).smoosh(['ph2','ph1','tau'],'transients'), alpha=0.2)
        fl.next('raw signal')
        fl.image(d)
        #fl.show();quit()
        d = d['t2':f_range]
        d = d['ph1',+1]['ph2',-2]
        d.ift('t2')
        fl.next('echoes')
        fl.plot(d.real,alpha=0.2,linewidth=0.5)
        fl.plot(abs(d),alpha=0.5,linewidth=1)
        #fl.show();quit()
        #fl.plot(1.76e-5, 'x', label='estimated signal at 0')
        NV = 250e-6*55.4*2*N_A # 400 μL, 55.4 M water molecs, 2 spins/molec
        nu0 = 14.89e6 # check me on this
        LambdaNMR = 1.55e-4 # 1 G/√W
        I = 0.5
        Vsignal = LambdaNMR * NV * (gammabar_H*2*pi) * I * (I+1) * (hbar*2*pi*nu0)**2 * sqrt(50)
        Vsignal /= 3 * k_B * (273+20)
        axhline(y=Vsignal,alpha=0.2)
        print("Vsignal expected",Vsignal)
        #fl.plot(Vsignal, 'x', label='theoretical signal at 0')
        fl.show();quit()
