"""Check integral error calculation -- simple version
======================================================

This is a simpler version of the check integral error calculation.
Rather than generating 100 fake inversion recovery curves, it generates only 100 fake scans that are identical.
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import integrate_limits, integral_w_errors

# sphinx_gallery_thumbnail_number = 1

init_logging(level="debug")
fl = figlist_var()
t2 = nddata(r_[0:1:1024j], "t2")
ph1 = nddata(r_[0, 2] / 4.0, "ph1")
ph2 = nddata(r_[0:4] / 4.0, "ph2")
signal_pathway = {"ph1": 0, "ph2": 1}
excluded_pathways = [(0, 0), (0, 3)]
# this generates fake clean_data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 10 Hz
clean_data = 21*exp(+1j*2*pi*100*t2 - t2*10*pi)
clean_data *= exp(signal_pathway["ph1"]*1j*2*pi*ph1)
clean_data *= exp(signal_pathway["ph2"]*1j*2*pi*ph2)
clean_data["t2":0] *= 0.5
fake_data_noise_std = 2.0
clean_data.reorder(["ph1", "ph2"])
bounds = (0, 200)  # seem reasonable to me
result = 0
n_repeats = 100 
all_manual_results = ndshape(clean_data) + (n_repeats, "repeats")
all_manual_results.pop("t2").pop("ph1").pop("ph2")
all_manual_results = all_manual_results.alloc().setaxis("repeats",r_[0:n_repeats]).set_error(0)
all_autointegrate_results = all_manual_results.C
print("shape of all results", ndshape(all_manual_results))
for j in range(n_repeats):
    data = clean_data.C
    data.add_noise(fake_data_noise_std)
    # at this point, the fake data has been generated
    data.ft(["ph1", "ph2"])
    # {{{ usually, we don't use a unitary FT -- this makes it unitary
    data /= 0.5 * 0.25  # the dt in the integral for both dims
    data /= sqrt(ndshape(data)["ph1"] * ndshape(data)["ph2"])  # normalization
    # }}}
    dt = diff(data.getaxis("t2")[r_[0, 1]]).item()
    data.ft("t2", shift=True)
    # {{{
    data /= sqrt(ndshape(data)["t2"]) * dt
    print(ndshape(data))
    error_pathway = (set(((j,k) for j in range(ndshape(data)['ph1']) for k in range(ndshape(data)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
    s_int,frq_slice = integral_w_errors(data,signal_pathway,error_pathway,
            fl=fl,return_frq_slice=True)
    # }}}
    manual_bounds = data["ph1", 0]["ph2", 1]["t2":frq_slice]
    N = ndshape(manual_bounds)["t2"]
    df = diff(data.getaxis("t2")[r_[0, 1]]).item()
    manual_bounds.integrate("t2")
    # N terms that have variance given by fake_data_noise_std**2 each multiplied by df
    print("#%d"%j)
    print(ndshape(data))
    std_off_pathway = (
        data["ph1", 0]["ph2", 0]["t2":bounds]
        .C.run(lambda x: abs(x)**2/2) # sqrt2 so variance is variance of real
        .mean_all_but(["t2"])
        .mean("t2")
        .run(sqrt)
    )
    print(std_off_pathway)
    propagated_variance_from_inactive = N * df ** 2 * std_off_pathway ** 2
    print(propagated_variance_from_inactive)
    quit()
    manual_bounds.set_error(array([sqrt(propagated_variance_from_inactive.data)]))
    all_manual_results["repeats", j] = manual_bounds
    all_autointegrate_results["repeats", j] = s_int
for thislabel, thisdata in [
        ("all autointegrate",all_autointegrate_results),
        ("all manual",all_manual_results),
        ]:
    fl.next(thislabel)
    fl.plot(thisdata,'o')
    thisavg = thisdata.data.mean()
    thisstd = std(thisdata.data.real)
    axhline(thisavg, color='k', ls=':', alpha=0.5)
    axhline(thisavg-thisstd, color='k', alpha=0.5)
    axhline(thisavg+thisstd, color='k', alpha=0.5)
fl.show()

