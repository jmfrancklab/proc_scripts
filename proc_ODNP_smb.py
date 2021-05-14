# -*- coding: utf-8 -*-
"""
ODNP Analysis
"""
# {{{ Imports
from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.correlation_alignment_ODNP import correl_align
from sympy import symbols, Symbol, latex,limit,init_printing
from matplotlib import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from sympy import exp as s_exp
# }}}
logger = init_logging("info")
fl = fl_mod()
t2 = symbols('t2')

def select_pathway(s,pathway):
    retval = s
    for k,v in pathway.items():
        retval = retval[k,v]
    return retval
    def as_scan_nbr(s):
        return s.C.setaxis('power','#').set_units('power','scan #')
                           
signal_pathway = {'ph1':1,'ph2':-2}
# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be


file_list = ['210422_210422_T177R1a_pR_KI_ONDP']
# For now, can only process one data set at a time because need to be able 
    #to save and load a list of T1 and power values to use the processing...
    # Manually enter for now:
T1_vals = r_[]
power_vals = []

# leave this as a loop, so you can load multiple files
for i,(searchstr, exp_type, nodename,
     postproc, freq_range, t_range, n_powers) in enumerate(zip(
        file_list,['Sam']*1,['signal']*1,
        ['spincore_ODNP_v1']*1,
        [(-400,400)]*1,
        [(None,0.03)]*1,
        [20]*1
        ):
    s = find_file(searchstr, exp_type = exp_type, expno = nodename,
                  postproc = postproc, lookup = postproc_dict, fl=fl)
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# In[]
    