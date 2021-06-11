"""
Check Integration
=================

Makes sure that automatically chosen integral bounds perform similar to or
better than what you would choose by hand.
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import integrate_limits, integral_w_errors
from numpy.random import seed

seed(2021)  # so the same random result is generated every time -- 2021 is meaningless
init_logging(level="debug")
fl = figlist_var()
t2 = nddata(r_[0:1:1024j], "t2")
vd = nddata(r_[0:1:40j], "vd")
ph1 = nddata(r_[0, 2] / 4.0, "ph1")
ph2 = nddata(r_[0:4] / 4.0, "ph2")
signal_pathway = {"ph1": 0, "ph2": 1}
excluded_pathways = [(0, 0), (0, 3)]
bounds = (60,140) # seem reasonable to me
    #(-90.91113281, 307.69921875),  # what's currently picked by the automatic routine

# this generates fake data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 10 Hz
data = 21*(1 - 2*exp(-vd / 0.2))*exp(+1j*2*pi*100*t2 - t2*10*pi)
data *= exp(signal_pathway["ph1"]*1j*2*pi*ph1)
data *= exp(signal_pathway["ph2"]*1j*2*pi*ph2)
data["t2":0] *= 0.5
fake_data_noise_std = 2.0
#data.add_noise(fake_data_noise_std)
data.reorder(["ph1", "ph2", "vd"])
# at this point, the fake data has been generated
data.add_noise(fake_data_noise_std)
data.ft(["ph1", "ph2"])
fl.next("what does a usual error bar look like?")
just_noise = nddata(r_[0:1:50j], "t")
just_noise.data *= 0
just_noise.add_noise(fake_data_noise_std)
just_noise.set_error(fake_data_noise_std)
fl.plot(just_noise, ".", capsize=6)
# {{{ usually, we don't use a unitary FT -- this makes it unitary
data /= 0.5 * 0.25  # the dt in the integral for both dims
data /= sqrt(ndshape(data)["ph1"] * ndshape(data)["ph2"])  # normalization
# }}}
dt = diff(data.getaxis('t2')[r_[0,1]]).item()
data.ft('t2',shift=True)
data /= sqrt(ndshape(data)['t2']) * dt
error_pathway = (set(((j,k) for j in range(ndshape(data)['ph1']) for k in range(ndshape(data)['ph2'])))
        - set(excluded_pathways)
        - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
s_int,frq_slice,mystd = integral_w_errors(data,signal_pathway,error_pathway,
        indirect='vd', fl=fl,return_frq_slice=True)
logger.debug(strm("check the std after FT", std(data["ph1", 0]["ph2", 0].data.real)))
# the sqrt on the next line accounts for the var(real)+var(imag)
fl.next("compare manual vs. automatic", legend=True)
# here I plot w/ manually chosen integration bounds:
manual_bounds = data["ph1", 0]["ph2", 1]["t2":frq_slice]
N = ndshape(manual_bounds)["t2"]
df = diff(data.getaxis("t2")[r_[0, 1]]).item()
logger.debug(
    strm(ndshape(manual_bounds), "df is", df, "N is", N, "N*df is", N * df)
)
manual_bounds.integrate("t2")
std_off_pathway = (
    data["ph1", 0]["ph2", 0]["t2":bounds]
    .C.run(lambda x: abs(x)**2/2)
    .mean_all_but(["t2","vd"])
    .mean("t2")
    .run(sqrt)
)
logger.debug(
    strm(
        "here is the std calculated from an off pathway",
        std_off_pathway,
        "does it match",
        fake_data_noise_std,
        "?",
    )
)
# N terms that have variance given by fake_data_noise_std**2 each multiplied by df
# the 2 has to do w/ real/imag/abs -- see check_integration_error
propagated_variance_from_inactive = N * df **2 *std_off_pathway ** 2
propagated_variance = N * df ** 2 * fake_data_noise_std ** 2
fl.plot(s_int,'.',capsize=6,label='std calculated from integral_w_error')
logger.debug(
    strm("manually calculated integral error is", sqrt(propagated_variance))
)
manual_bounds.set_error(sqrt(propagated_variance))
fl.plot(
    manual_bounds,
    ".",
    capsize=6,
    label=r"bounds $%4g\rightarrow%4g$" % bounds,
    alpha=0.5,
)
manual_bounds.set_error(sqrt(propagated_variance_from_inactive.data))
fl.plot(manual_bounds,
        '.',capsize=6,
        label = r'propagated from inactive std',alpha=0.5)
fl.plot(s_int,'.',capsize=6,label='error calculated with integral_w_errors',alpha=0.5)
fl.show()
