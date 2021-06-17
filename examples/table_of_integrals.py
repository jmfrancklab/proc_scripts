"""Convert 2D to Integrals with Errors
===================================

Take a 2D dataset and convert it to a table of integrals with errors, utilizing
all the bells and whistles (frequency and time selection, alignment, etc.)

Demonstrate on a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise, and with fluctuating field
(normally distributed field variation).
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import integrate_limits, integral_w_errors
from numpy.random import normal, seed
import sympy as s
from sympy import Symbol as sympy_symbol
from collections import OrderedDict
seed(2021)
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 4

# {{{ generate the fake data
init_logging(level="debug")
fl = figlist_var()
class fl_dummy_class (object):
    def plot(*args):
        pass
    def next(*args):
        pass
    def push_marker(*args):
        pass
    def pop_marker(*args):
        pass
fl_dummy = fl_dummy_class
def fake_data(
        expression,
        axis_coords,
        signal_pathway,
        direct = 't2',
        SD_sigma = [0.05,0.003],
        SD_amp = [1,5],
        scale=100,
        fake_data_noise_std = 1.0,
        fl=fl_dummy):
    """Generate fake data subject to noise and frequency variation.
    Parameters
    ==========
    expression: sympy expression
        Gives the functional form of the data.
    axis_coords: OrderedDict
        Gives nddata objects providing all the axis coordinates.
        **Very importantly**, these must be listed in the loop nesting order
        (outside in) in which they occur in the pulse program,
        or the frequency drift will not be modeled correctly.

        To enable *simulating echo-like data*, you can specify a direct axis
        that starts at a negative number.
        If you do this, the beginning of the axis will be re-set to 0 before returning.
    signal_pathway: dict
        Gives the signal pathway, with keys being phase cycling dimensions, and
        values being the corresponding Δp.
    scale: float (default 100) 
        amplitude of frequency variation
    """
    # this generates fake clean_data w/ a T₂ of 0.2s
    # amplitude of 21, just to pick a random amplitude
    # offset of 300 Hz, FWHM 10 Hz
    mysymbols = expression.atoms(sympy_symbol)
    missing = (
            set(str(j) for j in mysymbols)
            - set(axis_coords.keys()))
    assert len(missing) == 0, "all non-phase cycling symbols in your expression must have matching axis coordinates in the dictionary -- you are missing %s!"%str(missing)
    thefunction = lambdify(mysymbols, expression, 'numpy')
    clean_data = thefunction(*tuple(axis_coords[str(j)] for j in mysymbols))
    for j in signal_pathway.keys():
        clean_data *= exp(signal_pathway[j]*1j*2*pi*axis_coords[j])
    ## {{{ model frequency drift
    indirect_size = prod([ndshape(clean_data)[j] for j in axis_coords.keys() if j != direct])
    frq_noise = normal(scale=scale,size=indirect_size)
    frq_noise = frq_noise + 1j*normal(scale=scale,size=frq_noise.size)
    # {{{ control the spectral density of the shifts to be gaussian
    frq_noise = nddata(frq_noise,[-1],['temp'])
    N = ndshape(frq_noise)['temp']
    frq_noise.setaxis('temp',-0.5+r_[0:N]/N).set_units('temp','cycperscan')
    SD_gen = zip(SD_sigma,SD_amp)
    sigma,A = next(SD_gen)
    frq_noise_dens = A*exp(-frq_noise.fromaxis('temp')**2/2/sigma**2)
    for sigma,A in SD_gen:
        frq_noise_dens += A*exp(-frq_noise.fromaxis('temp')**2/2/sigma**2)
    frq_noise *= frq_noise_dens
    fl.push_marker()
    fl.next('frq-noise density')
    fl.plot(frq_noise)
    frq_noise.ift('temp')
    frq_noise /= sqrt(ndshape(frq_noise)['temp']) * frq_noise.get_ft_prop('temp','df') # normalization
    fl.next('frq-noise time domain')
    fl.plot(frq_noise)
    # }}}
    frq_noise = nddata(frq_noise.data.real,
            [ndshape(clean_data)[j] for j in axis_coords.keys() if j != direct],
            [j for j in axis_coords.keys() if j != direct])
    ## }}}
    data = clean_data.C
    data.add_noise(fake_data_noise_std)
    data *= exp(1j*2*pi*frq_noise*(data.fromaxis('t2'))) # the frequency shift
    # at this point, the fake data has been generated
    for j in signal_pathway.keys():
        data.ft(j)
    fl.pop_marker()
    data.setaxis(direct, lambda x: x-data.getaxis(direct)[0])
    data.ft(direct, shift=True)
    data.ift(direct)
    data.register_axis({direct:0})
    return data
t2, td, vd, ph1, ph2 = s.symbols('t2 td vd ph1 ph2')
echo_time = 5e-3
data = fake_data(
    21*(1 - 2*s.exp(-vd / 0.2))*s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi),
    OrderedDict([
        ("vd" , nddata(r_[0:1:40j], "vd")),
        ("ph1" , nddata(r_[0, 2] / 4.0, "ph1")),
        ("ph2" , nddata(r_[0:4] / 4.0, "ph2")),
        ("t2" , nddata(r_[0:0.2:256j]-echo_time, "t2"))]),
        {"ph1": 0, "ph2": 1})
data.reorder(["ph1", "ph2", "vd"])
fl.next("fake data -- time domain")
fl.image(data)
data.ft("t2")
# {{{ make it unitary again
data /= sqrt(ndshape(data)["t2"]) * data.get_ft_prop('t2','dt')
# }}}
fl.next("fake data -- freq domain")
fl.image(data)
# }}}
fl.show()
