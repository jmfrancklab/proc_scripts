from pyspecdata import *
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from scipy.optimize import minimize,curve_fit,least_squares
rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]

fl = figlist_var()
for searchstr,exp_type,expno,postproc,manual_phcyc,w0 in [
        #('ag_oct182019_bulk',1,True,0),
        #('ag_oct182019_w0_1',1,True,1),
        ('w8_200224','test_equip',1,'zg2h',
            True,3),
        #('ag_oct182019_w0_6',1,True,6),
        #('ag_oct182019_w0_8',1,True,8),
        #('ag_oct182019_w0_12',1,True,12),
        #('ag_bulk_d20',1,True,15),
        #('ab_jul092019_rtAG',2,False,8),
        #('jul012019_rtAG',2,False,12),
        #('ag_jul182019_w0_1',1,True,1),
        #('ag_jul182019_w0_3',1,True,3),
        #('ag_jul182019_w0_8',1,True,8),
        #('ag_jul312019_w0_8',3,True,8),#with paropt
        #('ag_jul312019_w0_12',1,True,12),#with paropt
        #('ag_aug072019_w0_1',1,True,1),
        #('ag_aug072019_w0_3_2h',3,True,3),
        #('ag_aug192019_w0_1_300_mM',1,True,1), #300 mM sample
        #('ag_aug192019_w0_3',1,True,3),
        #('ag_aug192019_w0_3_300_mM',1,True,3),#300 mM same water content as water loading of 12
        #('ag_aug232019_w0_1_900_mM',1,False,1),
        #('ag_sep062019_w0_3_IR',2,False,3)
        # these are the names of the data folders
        ]:
    label_id = '$w_0=%f$'%w0
    if manual_phcyc:
        s = find_file(searchstr,exp_type=exp_type,
                expno=expno, postproc=postproc,
                lookup=postproc_dict,dimname='ph') # here I am loading the data from the Google drive
        s.labels('ph',r_[0:4]).ft('ph') # assume it's a four step phase cycle
        s = s['ph',-1] # assuming it's a ninety pulse experiment
    else:
        s = find_file(searchstr, exp_type=exp_type, 
                expno=expno, postproc=postproc,
                lookup=postproc_dict,dimname='ph')
# here I am loading the data from the Google drive
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

