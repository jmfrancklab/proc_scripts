"""
Demonstrate Integrate Limits
============================

For this demonstration, we generate inversion
recovery data for a single peak, with a relatively
mild frequency variation, so that no serious
alignment is required before integration. We mimic
the 8-step phase cycle used for echo detection in
these experiments, and include the effect of the
echo time on the data detected in the time domain.

We use integrate_limits to detect the frequency
limits used for peak integration, based on a
matched Lorentzian filter on our frequency domain
data.

We illustrate the position of the frequency
limits with vertical lines on the final plot.

"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from numpy.random import normal, seed
from numpy.linalg import norm
import sympy as s
from collections import OrderedDict
seed(2021)
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 1
init_logging(level="debug")

with figlist_var() as fl:
    # {{{ generate the fake data
    # this generates fake clean_data w/ a T1 of 0.2s
    # amplitude of 21, just to pick a random amplitude
    # offset of 300 Hz, FWHM 10 Hz
    t2, td, vd, ph1, ph2 = s.symbols('t2 td vd ph1 ph2')
    echo_time = 5e-3
    data = fake_data(
        21*(1 - 2*s.exp(-vd / 0.2))*s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi),
        OrderedDict([
            ("vd" , nddata(r_[0:1:40j], "vd")),
            ("ph1" , nddata(r_[0, 2] / 4.0, "ph1")),
            ("ph2" , nddata(r_[0:4] / 4.0, "ph2")),
            ("t2" , nddata(r_[0:0.2:256j]-echo_time, "t2"))]),
            {"ph1": 0, "ph2": 1},
            scale=20.)
    # {{{ just have the data phase (not testing phasing here)
    data.setaxis('t2', lambda x: x-echo_time).register_axis({"t2":0})
    data = data['t2',0:-3] # dropping the last couple points avoids aliasing
    #                        effects from the axis registration
    #                        (otherwise, we get "droop" of the baseline)
    # }}}
    data.reorder(["ph1", "ph2", "vd"])
    fl.next("fake data -- time domain")
    fl.image(data)
    fl.next("FID sliced -- time domain")
    data = data['t2':(0,None)]
    data['t2',0] *= 0.5
    ph0 = data['t2',0].data.mean()
    ph0 /= abs(ph0)
    data /= ph0
    fl.image(data)
    data.ft("t2")
    fl.next("fake data -- freq domain")
    fl.image(data)
    for method in ['Lorentzian','Gaussian']:
        fl.basename = method + " filter:"
        freq_lim = integrate_limits(data['ph1',0]['ph2',1],
                convolve_method=method,
                fl=fl)
        fl.next("fake data -- show freq limit selection")
        fl.plot(data['ph1',0]['ph2',1])
        axvline(x=freq_lim[0])
        axvline(x=freq_lim[-1])
        print("Determined frequency limits via",method,"filter of",freq_lim)
    # }}}
    # }}}
