from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from hermitian_function_test import zeroth_order_ph
from scipy.optimize import minimize,curve_fit,least_squares
rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]

fl = figlist_var()
def calc_baseline(this_d,
        ph1lim,
        npts=5,
        guess=None,
        show_plots=True):
    if show_plots: fl.next('try baseline correction')
    if show_plots: fl.plots(this_d,
            label='before')
    this_d_tdom = this_d.C.ift('t2')
    blank_tdom = this_d_tdom.C
    blank_tdom.data[:] = 0
    def ve_to_params(ini_vec):
        phcorr0,phcorr1 = ini_vec[:2]
        baselin_vec = ini_vec[2:].view(complex128)
        return phcorr0,phcorr1,baseline_vec
    def apply_corr(ini_vec):
        phcorr0,phcorr1,baseline_vec = vec_to_params(ini_vec)
        d_test = this_d.C
        d_test *= exp(-1j*phcorr1*d.fromaxis('t2')-1j*phcorr0)
        d_test.ift('t2')
        retval = d_test['t2',0:len(baseline_vec)] = retval
        return d_test.ft('t2')
    def generat_baseline(baseline_vec):
        baseline_data = blank_tdom.C
        baseline_data['t2',0:len(baseline_vec)] += baseline_vec
        return baseline_data
    def costfun(ini_vec):
        d_test = apply_corr(ini_vec)
        return abs(d_test.real).sum('t2').data.item()
    max_val = abs(this_d_tdom.data).max()
    print(max_val)
    mybounds = r_[-max_val,max_val][newaxis,:]*ones(npts*2)[:,newaxis]
    print(shape(mybounds))
    print(mybounds)
    mybounds = r_[
            r_[-pi,pi,-ph1lim,ph1lim].reshape(-1,2),
            mybounds]
    if guess is None:
        guess = zeros(npts*2+2)
    else:
        guess = r_[guess[0].real,
                guess[1].real,
                guess[2:].view(float64)]
        res = minimize(costfun, (guess,),
                method='L-BFGS-B',
                bounds=mybounds,
                )
        phcorr0,phcorr1,baseline_vec = vec_to_params(res.x)
        baseline = generate_baseline(baseline_vec)
        if show_plots: fl.plot(this_d*exp(-1j*phcorr1*d.fromaxis('t2')-1j*phcorr0)+baseline.C.ft('t2'))
        return phcorr0,phcorr1,baseline

for exp_name,expno,manual_phcyc,w0 in [
        #('ag_oct182019_bulk',1,True,0),
        #('ag_oct182019_w0_1',1,True,1),
        ('ag_oct182019_w0_3',1,False,3),
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
    print((exp_name,expno))
    label_id = '$w_0=%f$'%w0
    if manual_phcyc:
        s = find_file(exp_name, exp_type='NMR_Data_AG', dimname='ph', expno=expno) # here I am loading the data from the Google drive
        s.labels('ph',r_[0:4]).ft('ph') # assume it's a four step phase cycle
        s = s['ph',-1] # assuming it's a ninety pulse experiment
    else:
        s = find_file(exp_name, exp_type='NMR_Data_AG', expno=expno) # here I am loading the data from the Google drive
    dw = diff(s.getaxis('t2')[0:2]).item()
    # {{{ determine the phase corrections
    def phasecorrect(s):
        fl.push_marker() 
        ph1 = nddata(r_[-5:5:70j]*dw,'phcorr')
        dx = diff(ph1.getaxis('phcorr')[r_[0,1]]).item()
        ph1 = exp(-1j*2*pi*ph1*s.fromaxis('t2'))
        s_cost = s * ph1
        ph0 = s_cost.C.sum('t2')
        ph0 /= abs(ph0)
        s_cost /= ph0
        fl.next('phasing cost function')
        s_cost.run(real).run(abs).sum('t2')
        fl.plot(s_cost,'.')
        ph1_opt = s_cost.argmin('phcorr').item()
        print('optimal phase correction',repr(ph1_opt))
        # }}}
        # {{{ apply the phase corrections
        def applyphase(arg,ph1):
            arg *= exp(-1j*2*pi*ph1*arg.fromaxis('t2'))
            ph0 = arg.C.sum('t2')
            ph0 /= abs(ph0)
            arg /= ph0
            return arg
        def costfun(ph1):
            if type(ph1) is ndarray:
                ph1 = ph1.item()
            temp = s.C
            retval = applyphase(temp,ph1).run(real).run(abs).sum('t2').item()
            return retval
        print("rough opt cost function is",costfun(ph1_opt))
        r = minimize(costfun,ph1_opt,
                bounds=[(ph1_opt-dx,ph1_opt+dx)])
        assert r.success
        s = applyphase(s,r.x.item())
        fl.plot(r.x,r.fun,'x')
        fl.pop_marker()
        return s
    s.ft('t2', shift=True)
    s = phasecorrect(s)
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

