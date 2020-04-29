from pyspecdata import *
from Utility import dBm2power
#to use type s = load_data("nameoffile",exptype,expno)
def load_data_bruker(searchstr,exptype,expno):
    if (exptype=='IR_noecho'):
        files = search_filename(searchstr, 'test_equip')
        s = find_file(searchstr, exp_type='test_equip', dimname = 'indirect', expno=expno)
        fl = figlist_var()
        print(ndshape(s))
        print(s.getaxis('indirect'))
        s.chunk('indirect',['indirect','ph1','ph2'],[-1,4,2]) #expands the indirect dimension into indirect, ph1, and ph2. inner most dimension is the inner most in the loop in pulse sequence, is the one on the farthest right. brackets with numbers are the number of phase cycle steps in each one. the number of steps is unknown in 'indirect' and is therefore -1.
        print(s.getaxis('indirect'))
        print(s.getaxis('ph1'))
        print(s.getaxis('ph2'))
        s.setaxis('ph1',r_[0:4.]/4) #setting values of axis ph1 to line up
        s.setaxis('ph2',r_[0:2.]/4) #setting values of axis ph1 to line up
        s.setaxis('indirect', s.get_prop('vd'))
        fl.next('phased coherence domain') #switch to time domain as a string based name for fig
        s.ft(['ph1','ph2']) #fourier transforming from phase cycle dim to coherence dimension
        s.reorder(['indirect','t2'], first=False)
        fl.image(s) #labeling
        return s

