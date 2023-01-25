import numpy as np
from numpy import r_
from pyspecdata import *
from pyspecProcScripts import *
"""Here, we want to calculate the time-domain variance to use in error propagation.
But, to make sure we calculate only noise, we want to mask out portions of the frequency 
domain.
We propose that if we use a unitary transform,
Parseval's theorem
tells us we can calculate :math:`\sigma_t^2` (time domain variance)
directly from the frequency-domain variance, *i.e.* :math:`\sigma_\nu^2=\sigma_t^2`.
To confirm this, we construct a "spectrum" of pure noise
and generate a frequency-masked noise,
and show that is the same as the 
unmasked time-domain.
"""
#{{{generate data with just noise with a phase cycling dimension and repeats dimension
N = 1024
n_repeats = 50
signal_window = (-100,200) # wherever my "peak" shows up
signal_pathway = {'ph':1}
example_data = nddata(np.random.normal(size=4*n_repeats*N)
        +1j*np.random.normal(size=4*n_repeats*N), [4,n_repeats,N], ['ph','repeats','t'])
example_data.setaxis('ph',r_[0:4]/4)
example_data.setaxis('repeats',r_[0:n_repeats])
example_data.setaxis('t',r_[0:1:1j*N])
#}}}
#{{{Test if data is unitary
# this data only has noise -- it should have a norm of 1
# the /2 comes from the variance along both real and imaginary
print("The original std of the data in the time domain is:",
        sqrt(example_data.C.run(np.var,'t').mean('ph')/2))
# the way that we do FT is parseval preserved?
temp = example_data.C.ft('t', shift=True)
print("If we apply just FT as we normally would the std in the frequency domain is:",
        sqrt(temp.run(np.var,'t').mean('ph')/2))
# it's not!  I need to use a unitary FT for this to work
print("These values are NOT the same so we need a unitary FT for this to work")
example_data.ft('t', shift=True, unitary=True)
example_data.ft('ph', unitary=True)
print("When we apply a unitary FT the std of the frequency is:",
        sqrt(example_data.C.run(np.var,'t').mean('ph')/2))
# now, I can just calculate the "time domain" noise variance in the
# frequency domain, where it's easier to mask out regions of the coherence
# domain where I expect there is signal (or phase cycling noise)
#}}}
# {{{ I'm doing a mildly odd thing where I'm using "nan" to identify signal I
# want to exclude from the variance calculation -- i.e. to mask it.  This is
# assuming that I have signal that I'm not interested in including in the
# calculation.
temp = select_pathway(example_data, signal_pathway)
temp.data[:] = nan # note how I am NOT acting on a copy -- I am trying to
#                    manipulate the data at its original memory position!
# for the most complicated case I'll also say I want to exclude phase cycling
# noise -- so also exclude everything from the signal bandwidth
# this will give a conservative (small) estimate of the noise
temp = example_data['t':signal_window]
temp.data[:] = nan
with figlist_var(black=True) as fl:
    fl.next('show the mask in white')
    forplot = example_data.C
    # in pyspecdata, nan shows up as the opposite (black vs. white) color vs. 0
    fl.image(forplot)
#}}}    
#{{{ Calculate the std using new functions
# now, I can do this:
example_data.run(masked_var_multi, 't')
example_data.run(masked_mean_multi,'ph')
# It only includes datapoints I have not
# "masked", so it works. Note that the following result is the time-domain
# variance that does *not* include regions of the coherence domain where we
# expect to find signal (or the effect of phase cycling noise)
print("The std when using the mask on unitary data is:",
        sqrt(example_data.data))
#}}}
