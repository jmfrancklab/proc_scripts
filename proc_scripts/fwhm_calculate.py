from pylab import subplots
import numpy as np
import logging
def fwhm_calculator(self, axis='t2',fl=None):
    signal_sign = self.C.sum(axis).run(np.real).run(np.sign)
    temp = s.real * signal_sign
    # pulled from apodization code
    sigma = nddata(np.linspace(1e-5,1e3,1000),'sigma').set_units('sigma','s')
    s_avg = self.C.mean_all_but('t2')
    gaussians = np.exp(-s_avg.C.fromaxis('t2')**2/2/sigma**2)
    signal_E = (abs(s_avg * gaussians)**2).sum('t2')
    signal_E /= signal_E.data.max()
    filter_width = abs(signal_E-1/sqrt(2)).argmin('sigma').item()
    if fl is not None:
        fl.next('signal Energy')
        fl.plot(signal_E, human_units=False)
        fl.plot(signal_E['sigma':(filter_width,filter_width+1e-6)],'o', human_units=False)
    fwhm = filter_width
    fl.push_marker()
    if fl is not None:
        fig, (ax1,ax2) = subplots(2,1)
        fl.next("integration diagnostic", fig=fig)
        fl.plot(temp, ax=ax1)
    temp.mean_all_but(axis)
    if fl is not None:
        fl.plot(temp/abs(temp.data).max(), ax=ax2)
    return fwhm    

