"""
Phasing and Timing Correction
=============================

Take real data with varying echo times, 
and demonstrate how we can automatically find the zeroth order phase and the
center of the echo in order to get data that's purely real in the frequency
domain.

This example provides a challenging test case, with low SNR data (from AOT RMs),
one of which has a very short echo time.
"""
from pyspecdata import *
from pyspecProcScripts import *
from pylab import *
import sympy as s
from collections import OrderedDict
from numpy.random import normal
from scipy.signal import tukey
init_logging(level='debug')

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery

# sphinx_gallery_thumbnail_number = 1
t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
f_range = (-500, 500)
filename = '210604_50mM_4AT_AOT_w11_cap_probe_echo'
signal_pathway = {'ph1':1,'ph2':0}
with figlist_var() as fl:
    for nodename,file_location,postproc,label,alias_slop in [
        ('tau_1000','ODNP_NMR_comp/Echoes','spincore_echo_v1',
            'tau is 1 ms',1),
        ('tau_3500','ODNP_NMR_comp/Echoes','spincore_echo_v1',
            'tau is 3.5 ms',3),
        ('tau_11135','ODNP_NMR_comp/Echoes','spincore_echo_v1',
            'tau is 11.135 ms',3),
            ]:
        data = find_file(filename,exp_type=file_location,expno=nodename,
                postproc=postproc,lookup=lookup_table)
        fl.basename = "(%s)" % label
        fig, ax_list = subplots(1, 3, figsize = (7,7))
        fig.suptitle(fl.basename)
        data.reorder(['ph1','ph2','nScans','t2'])
        fl.next("Data processing", fig=fig)
        fl.image(data['t2':(-1e3,1e3)], ax=ax_list[0])
        ax_list[0].set_title("Raw Data")
        data = data['t2':f_range]
        data.ift("t2")
        best_shift = hermitian_function_test(
            select_pathway(data.C.mean("nScans"), signal_pathway),
            aliasing_slop=alias_slop,
            fl=fl
        )
        logging.info(strm("best shift is:",best_shift))
        data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2": 0})
        data /= zeroth_order_ph(select_pathway(data,signal_pathway))
        data.ft('t2')
        fl.image(data['t2':(-1e3,1e3)], ax=ax_list[1])
        ax_list[1].set_title("Phased and centered (ν)")
        data.ift("t2")
        fl.image(data, ax=ax_list[2], human_units=False)
        ax_list[2].set_title("Phased and Centered (t)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
