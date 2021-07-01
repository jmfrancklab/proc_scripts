"""Align data with significant frequency drift
==============================================

Takes a 2D data set and applies proper phasing corrections followed by 
aligning the data through a correlation routine.
"""
from pyspecdata import *
from pyspecProcScripts import *
from pylab import *
import sympy as s
from collections import OrderedDict
from numpy.random import normal,seed
seed(2021)
rcParams['image.aspect'] = 'auto' #needed for sphinx gallery

t2, td, vd, power, ph1, ph2 = s.symbols('t2 td vd power ph1 ph2')
echo_time = 1e-3
datasets = []
f_ranges = []
s_pathways = []
fl = figlist_var()

def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval    

#{{{generate IR fake data
#this generates fake clean data with a T2 of 0.2s
# amplitude of 23 (arbitrary), offset of 300 Hz,
# FWHM 10 Hz
def make_IR_data():
    IR_data= fake_data(
        23*(1 - 2*s.exp(-vd / 0.2)) * s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi),
        OrderedDict([
            ('vd', nddata(r_[0:1:40j],'vd')),
            ('ph1', nddata(r_[0:4]/4.0,'ph1')),
            ('ph2', nddata(r_[0,2]/4.0,'ph2')),
            ('t2', nddata(r_[0:0.2:256j]-echo_time,'t2'))]),
            {'ph1':0,'ph2':1})
    IR_data.reorder(["ph1", "ph2", "vd"])
    fl.next("IR Data -- Time Domain")
    fl.image(IR_data)
    f_range = (-400,400)
    signal_pathway = {'ph1':0,'ph2':1}
    IR_data.ft('t2')
    IR_data /= sqrt(ndshape(IR_data)['t2'])*IR_data.get_ft_prop('t2','dt')
    fl.next('IR Data -- Frequency Domain')
    fl.image(IR_data)
    return IR_data, f_range, signal_pathway
#}}}

#{{{generate fake enhancement data
def make_enhancement_data():
    C_SL = 150e-6
    E_f_range = (-200,400)
    E_signal_pathway = {'ph1':1}
    enhancement_data = fake_data(
            23*(1-(32*power/(0.25 +power))*C_SL*659.33)*s.exp(+1j*2*s.pi*100*(t2) - abs(t2)*50*s.pi) ,
            OrderedDict([
                ('power', nddata(r_[0:4:25j],'power')),
                ('ph1', nddata(r_[0:4]/4.0,'ph1')),
                ('t2', nddata(r_[0:0.2:256j]-echo_time,'t2'))]),
            {'ph1':1})
    enhancement_data.reorder(['ph1','power','t2'])
    fl.next('Enhancement Data -- Time Domain')
    fl.image(enhancement_data)
    enhancement_data.ft('t2')
    enhancement_data /= sqrt(ndshape(enhancement_data)['t2'])*enhancement_data.get_ft_prop('t2','dt')
    fl.next('Enhancement Data -- Frequency Domain')
    fl.image(enhancement_data)
    return enhancement_data,E_f_range,E_signal_pathway
    #}}}

#{{{Function that will apply phase corrections and the correlation alignment
def processing(data, f_range, signal_path,label='',indirect=None):
    #{{{Changing sign of all signal to be the same so we don't have to
    # deal with the null of the signal
    myslice = data['t2':f_range]
    mysgn = select_pathway(myslice,signal_path).real.sum('t2').run(np.sign)
    data *= mysgn
    #}}}
    data = data['t2':f_range]
    data.ift('t2')
    #{{{Applying the phase corrections
    best_shift,max_shift = hermitian_function_test(select_pathway(data,signal_path).C.mean(indirect))
    data.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('%s After hermitian function test -- Time domain'%label)
    fl.image(data)
    data.ft('t2')
    fl.next('%s After hermitian function test -- Frequency domain'%label)
    fl.image(data)
    data.ift('t2')
    ph0 = select_pathway(data,signal_path)['t2':0]
    ph0 /= abs(ph0)
    data /= ph0
    fl.next('%s After phasing corrections applied'%label)
    fl.image(data)
    #}}}
    #{{{Applying Correlation Routine to Align Data
    data.ft('t2')
    data,opt_shift,sigma = correl_align(data,indirect_dim=indirect,
            signal_pathway = signal_path,
            sigma = 50)
    data.ift('t2')
    for k,v in signal_path.items():
        data.ft([k])
    fl.next('%s Aligned Data -- Time Domain'%label)
    fl.image(data)
    data.ft('t2')
    fl.next('%s Aligned Data -- Frequency Domain'%label)
    fl.image(data)
    #}}}

#{{{Process the fake IR data and the fake enhancement data
Ep, f_range, signal_path = make_enhancement_data()
Ep = processing(Ep,f_range,signal_path,label='Enhancement',indirect='power')
IR,IR_f_range, IR_signal_path = make_IR_data()
IR = processing(IR, IR_f_range, IR_signal_path, label='IR',indirect='vd')
#}}}
fl.show()


