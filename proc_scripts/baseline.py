def calc_baseline(this_d,
        ph1lim, 
        npts=5,
        guess=None,
        show_plots=True):
    """
    
    """
    if show_plots: fl.next('try baseline correction')
    if show_plots:fl.plot(this_d,
            label='before')
    this_d_tdom = this_d.C.ift('t2')
    blank_tdom = this_d_tdom.C
    blank_tdom.data[:] = 0
    def vec_to_params(ini_vec):
        phcorr0,phcorr1 = ini_vec[:2]
        baseline_vec = ini_vec[2:].view(complex128)
        return phcorr0,phcorr1,baseline_vec
    def apply_corr(ini_vec):
        phcorr0,phcorr1,baseline_vec = vec_to_params(ini_vec)
        d_test = this_d.C
        d_test *= exp(-1j*phcorr1*d.fromaxis('t2')-1j*phcorr0)
        d_test.ift('t2')
        retval = d_test['t2',0:len(baseline_vec)] + baseline_vec
        d_test['t2',0:len(baseline_vec)] = retval
        return d_test.ft('t2')
    def generate_baseline(baseline_vec):
        baseline_data = blank_tdom.C
        baseline_data['t2',0:len(baseline_vec)] += baseline_vec
        return baseline_data
    def costfun(ini_vec):
        d_test = apply_corr(ini_vec)
        return abs(d_test.real).sum('t2').data.item()
    max_val = abs(this_d_tdom.data).max()
    print(max_val)
    mybound = r_[-max_val,max_val][newaxis,:]*ones(npts*2)[:,newaxis]
    print(shape(mybounds))
    print(mybounds)
    mybound = r_[
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
    if show_plots: fl.plot(this_d*exp(-1j*phcorr1*d.fromaxis('t2')-1j*phcorr0)+baseline.C.ft('t2'),
            label='after')
    return phcorr0,phcorr1,baseline


