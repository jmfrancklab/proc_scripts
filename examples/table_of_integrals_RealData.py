"""Convert 2D  Real Data into Integrals with Errors
===================================================

Take a 2D dataset and convert it to a table of integrals with errors, utilizing
all the bells and whistles (frequency and time selection, alignment, etc.)

Demonstrate on a fake dataset of an inversion recovery with multiple repeats (φ
× t2 × vd × repeats) w/ normally distributed random noise, and with fluctuating field
(normally distributed field variation).
"""
from pylab import *
from pyspecdata import *
from pyspecProcScripts import *
from numpy.random import normal, seed
import sympy as s
from collections import OrderedDict
rcParams['image.aspect'] = 'auto' # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 8

init_logging(level="debug")
fl = fl_mod()
t2, td, vd, power, ph1, ph2 = s.symbols('t2 td vd power ph1 ph2')
signal_pathway = {'ph1':1,'ph2':0}
t_range= (0,85e-3)

with figlist_var() as fl:
    for filename, exp_type, nodename, postproc, indirect, clock_correction, label, f_range in [
            ('201113_TEMPOL_capillary_probe_DNP_1','ODNP_NMR_comp/old', 'signal',
                'spincore_ODNP_v2','power',False,'no power', 
                (-0.5e3,0.5e3))
            ]:
        fl.basename = "(%s)"%label
        s = find_file(filename, exp_type=exp_type, expno=nodename,
                postproc=postproc,lookup=lookup_table,fl=fl)
        if 'nScans' in s.dimlabels:
            s.mean('nScans')
        myslice = s['t2':f_range]
        mysgn = determine_sign(select_pathway(myslice,signal_pathway))
        s_int, s = process_data(s,signal_pathway=signal_pathway,
                searchstr = label, f_range=f_range, t_range=t_range,
                sgn = mysgn,indirect='power', Real =True, alias_slop=3,
                clock_correction=clock_correction,
                fl=fl)
