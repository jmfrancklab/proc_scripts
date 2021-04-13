from pylab import *
from pyspecdata import *
from proc_scripts import integrate_limits, integral_w_errors
init_logging(level='debug')
fl=figlist_var()
t2 = nddata(r_[0:1:1024j], 't2')
vd = nddata(r_[0:1:15j], 'vd')
ph1 = nddata(r_[0,2]/4.,'ph1')
ph2 = nddata(r_[0:4]/4.,'ph2')
signal_pathway = {'ph1':0,'ph2':1}
excluded_pathways = [(0,0),(0,3)]
# this generates fake data w/ a T₂ of 0.2s
# amplitude of 21, just to pick a random amplitude
# offset of 300 Hz, FWHM 10 Hz
data = 21*(1-2*exp(-vd/0.2)) * exp(+1j*2*pi*100*t2-t2*10*pi)
data *= exp(signal_pathway['ph1']*1j*2*pi*ph1)
data *= exp(signal_pathway['ph2']*1j*2*pi*ph2)
data['t2':0] *= 0.5
fake_data_noise_std = 0.1
data.add_noise(fake_data_noise_std)
data.reorder(['ph1','ph2','vd'])
# at this point, the fake data has been generated
data.ft(['ph1','ph2'])
# {{{ usually, we don't use a unitary FT -- this makes it unitary
data /= 0.5*0.25 # the dt in the integral for both dims
data /= sqrt(ndshape(data)['ph1']*ndshape(data)['ph2']) # normalization
# }}}
print("check the std after FT",std(data['ph1',0]['ph2',0].data.real))
data.ft('t2', shift=True)
fl.next('before integration')
fl.image(data, alpha=0.5)
fl.next('manual limits', legend=True)
# here I plot w/ manually chosen integration bounds:
for bounds in [(0,200), # seem reasonable to me
    (-95.90625,295.7109375) # what's currently picked by the automatic routine
    ]:
    manual_bounds = data['ph1',0]['ph2',1]['t2':bounds]
    N = ndshape(manual_bounds)['t2']
    df = diff(data.getaxis('t2')[r_[0,1]]).item()
    print(ndshape(manual_bounds),"df is",df,"N is",N,"N*df is",N*df)
    manual_bounds.integrate('t2')
    # N terms that have variance given by fake_data_noise_std**2 each multiplied by df
    propagated_variance = N * df**2 * fake_data_noise_std**2
    manual_bounds.set_error(sqrt(propagated_variance))
    fl.plot(manual_bounds, '.', capsize=6,
            label=r'bounds $%4g\rightarrow%4g$'%bounds,
            alpha=0.5)
# Now you need to run your code that automatically chooses integration bounds
# and also assigns error
fl.next('automatic limits and error')
error_pathway = (set(((j,k) for j in range(2) for k in range(4)))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
error_pathway = [{'ph1':j,'ph2':k} for j,k in error_pathway]
data = integral_w_errors(data,signal_pathway,error_pathway)
fl.plot(data,'.',label='real', capsize=6)
fl.plot(data.imag,'.',label='imaginary', capsize=6)
fl.show()
