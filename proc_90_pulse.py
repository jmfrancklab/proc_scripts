from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from proc_scripts import zeroth_order_ph,calc_baseline,ph1_real_Abs
from proc_scripts.load_data import postproc_dict
from scipy.optimize import minimize,curve_fit,least_squares
rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]

fl=figlist_var()
for searchstr, exp_type, which_exp, postproc, manual_phcyc, w0 in [
        ('ag_oct182019_w0_3', 'test_equip', 1, 'zg2h', False, 3),
        ]:
    label_id='$w_0=%f$'w%0
    if manual_phcyc:
        s = find_file(searchstr, exp_type=exp_type, expno= which_exp,
                postproc=postproc,lookup=postproc_dict,
                dimname='ph')
        s.labels('ph',r_[0:4]).ft('ph') #assuming it is a four step cycle
        s = s['ph',-1] #assuming it is a 90 pulse exp
    else:
        s = find_file(searchstr, exp_type=exp_type, expno=which_exp,
                postproc=postproc, lookup=postproc_dict)
    dw = diff(s.getaxis('t2')[0:2]).item()
    # {{{ determine the phase corrections
    s.ft('t2', shift=True)
    s = ph1_real_Abs(s)
    # }}}
    s_conv = s.C
    s_conv.ift('t2')
    s_conv *= exp(-8*pi*s_conv.fromaxis('t2'))
    s_conv.ft('t2')
    fl.next('raw - frequency domain')
    fl.image(s)
    fl.next('raw - time domain')
    d.ift('t2')
    fl.image(s)
    fl.show();quit()
    fl.next('raw data-frequency domain')

    peak = abs(d_conv).contiguous(lambda x:
            x > 0.1*x.data.max())
    def filter_range(thisrange):
        mask = diff(thisrange, axis=1) > 0.1e-6*ones((1,2))
        thisrange = thisrange[mask].reshape((-1,2))
        return thisrange
    peak = filter_range(peak)
    assert peak.shape[0] == 1, "there should only be one peak here"+repr(peak.shape)
    peak = peak[0,:]
    fl.plot(abs(d_conv), alpha=0.1, color='b',
            linewidth=1)
    fl.plot(abs(d_conv['t2':tuple(peak)]), alpha=0.1, color='k',
            linewidth=2)
    peak_start, peak_end = peak
    # {{{ manually (using numpy, not pyspecdata) slice out everything but the peak
    # need to upgrade pyspecdata to allow something like
    # d_conv = d_conv['t2':[(None,peak_start),(peak_end,None)]]
    d_baseline = d.C
    xaxis = d_baseline.getaxis('t2')
    mask = logical_or(xaxis<peak_start, xaxis>peak_end)
    d_baseline.data = d_baseline.data[mask]
    d_baseline.setaxis('t2',xaxis[mask])
    # }}}

    fl.next('region for baseline', legend=True)
    fl.plot(d_baseline.real,'.', alpha=0.5, label='real')
    fl.plot(d_baseline.imag,'.', alpha=0.5, label='imag')
    grab_y = gca().get_ylim()
    fl.plot(d.real, alpha=0.1)
    fl.plot(d.imag, alpha=0.1)
    c,_ = d_baseline.real.polyfit('t2', order=5)
    baseline = d.fromaxis('t2').run(lambda x: sum(
        c[j]*x**j for j in range(len(c))))
    fl.plot(baseline, alpha=0.5, linewidth=5, label='baseline')
    ylim(grab_y)
    fl.next('time domain of baseline correction (for jmf info)')
    d.ift('t2')
    fl.plot(d, alpha=0.5, label='phasing only')
    d.ft('t2')
    fl.next('baseline correction', legend=True)
    fl.plot(d, alpha=0.5, label='phasing only')
    d -= baseline
    fl.plot(d, alpha=0.5, label='phasing + baseline')
    fl.next('time domain of baseline correction (for jmf info)')
    d.ift('t2')
    fl.plot(d, alpha=0.5, label='phasing + baseline')
    d.ft('t2')
    fl.next('baseline correction')
    d = phasecorrect(d)
    fl.plot(d, alpha=0.5, label='phasing + baseline + phasing')
    ylim(grab_y)
    fl.next('time domain of baseline correction (for jmf info)')
    d.ift('t2')
    fl.plot(d, alpha=0.5, label='phasing + baseline + phasing')
    d.ft('t2')
    # now, you want to try to redo the phase correction procedure (wrap it as a function, so you can reuse it)

fl.show() # showing the plot
quit()

    
