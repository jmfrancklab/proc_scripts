from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from proc_scripts import zeroth_order_ph, calc_baseline, ph1_real_Abs
from proc_scripts.load_data import postproc_dict
from scipy.optimize import minimize,curve_fit,least_squares
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]
fl = figlist_var()
for searchstr, exp_type, which_exp, postproc, manual_phcyc, fl in [
        ('w8_200917','test_equip',5,'ag_zg2h',True,fl),
        ]:
    s = find_file(searchstr, exp_type=exp_type,
                expno=which_exp, postproc=postproc, lookup=postproc_dict,fl=fl)  
    #fl.show();quit()
    if manual_phcyc:
        fl.basename = searchstr
        #fl.next('image of phase cycle domain')
        #fl.image(s)
        print(ndshape(s))
        fl.show();quit()
        s.setaxis('ph',r_[0,1,2,3]/4.)
        s.ft('ph')# based on 90 pulse experiment, we seem to want ft rather than ift for deuterium (phcyc goes around the circle the wrong way)
        fl.next('image of coherence domain')
        fl.image(s)
        s = s['ph',-1].C
    fl.show();quit()
    s.ft('t2')
    fl.next('before phase correction')
    fl.plot(s)
    fl.next('after phase correction')
    dw = diff(s.getaxis('t2')[r_[0,-1]]).item()
    s = ph1_real_Abs(s,dw,fl=fl)
    s_before = s.C
    fl.plot(s)
    #fl.show();quit()
    #{{{slice and apply polynomial method of baseline
    peak = abs(s).contiguous(lambda x:
            x > 0.05*x.data.max())
    def filter_range(thisrange):
        mask = diff(thisrange, axis=1) > 0.05e-6*ones((1,2))
        thisrange = thisrange[mask].reshape((-1,2))
        return thisrange
    peak = filter_range(peak)
    #assert peak.shape[0] == 1, "there should only be one peak here"+repr(peak.shape)
    peak = peak[0,:]
    peak_start,peak_end = peak
    baseline = s.C
    xaxis = baseline.getaxis('t2')
    mask = logical_or(xaxis<peak_start, xaxis>peak_end)
    baseline.data = baseline.data[mask]
    baseline.setaxis('t2',xaxis[mask])
    fl.next('region for baseline',legend=True)
    fl.plot(baseline.real,'.',alpha=0.5, label='real')
    fl.plot(baseline.imag,'.',alpha=0.5,label='imag')
    grab_y = gca().get_ylim()
    fl.plot(s.real,alpha=0.1)
    fl.plot(s.imag, alpha=0.1)
    c,_ = baseline.real.polyfit('t2',order=5)
    baseline = s.fromaxis('t2').run(lambda x: sum(
        c[j]*x**j for j in range(len(c))))
    #}}}
    #{{{baseline function application
    #phcorr0,phcorr1,baseline = calc_baseline(this_d=s,ph1lim=50,npts=5,guess=None,fl=fl)
    #fl.next('baseline to be applied',legend=True)
    #fl.plot(baseline.real,alpha=0.5,label='baseline real')
    #fl.plot(baseline.imag,alpha=0.5,label='baseline imaginary')
    #fl.plot(s, alpha=0.1,label='actual data without baseline')
    #s = baseline
    #}}}

    fl.next('baseline subtracted',legend=True)
    fl.plot(s,alpha=0.5,label='after baseline application')
    fl.plot(s_before,alpha=0.5,label='before baseline is applied')
    fl.show();quit()
    # }}}
    s = ph1_real_Abs(s,dw,fl=fl)
    fl.next('with additional first order phasing',legend=True)
    fl.plot(baseline)
    fl.plot(s_before)


    # {{{ normalize, shift the axis, and select the center 300 Hz of the spectrum
    fl.next('normalized and centered')
    max_val_at = s.C.argmax('t2').item() #adding maxes to be put on git screen
    s /= s['t2':max_val_at]
    s.setaxis('t2',lambda x: x-max_val_at)
    print("SFO1",s.get_prop('acq')['SFO2'])
    fl.plot(s['t2':(-150,150)], alpha=0.5, label=r"%s $\nu_0$=%f ppm"%(searchstr,max_val_at/s.get_prop('acq')['SFO1']),human_units=False)
    # }}}
fl.show() # showing the plot
quit()
