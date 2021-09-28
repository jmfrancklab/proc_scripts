"""
Phasing and Timing Correction
=============================

Take real data with varying echo times, 
and demonstrate how we can automatically find the zeroth order phase and the
center of the echo in order to get data that's purely real in the frequency
domain.

Here, we specifically check to see whether or not the offset between the
programmed τ and the center of the echo
(as found by :func:`hermitian_function_test`)
is consistent.
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
f_range = (-200, 200)
filename = '201113_TEMPOL_capillary_probe_var_tau_1'
signal_pathway = {'ph1':1,'ph2':0}
est_tau = []
diff_tau = []
with figlist_var() as fl:
    for nodename,file_location,postproc,label in [
        ('var_tau','ODNP_NMR_comp/var_tau','spincore_var_tau_v1',
            'tau is 1 ms'),
            ]:
        data = find_file(filename,exp_type=file_location,expno=nodename,
                postproc=postproc,lookup=lookup_table)
        tau_list = list(data.getaxis('tau'))
        data.reorder(['ph1','ph2','tau','t2'])
        data = data['t2':f_range]
        for j in range(len(tau_list)):
            programmed_tau = tau_list[j]
            logger.info(strm("programmed tau:",programmed_tau))
            tau_data = data['tau',j]
            this_data = tau_data.C
            this_data.ift("t2")
            if programmed_tau > 0.04:
                print("extending data")
                this_data.extend("t2", tau_list[-1]*2)
            best_shift = hermitian_function_test(
                select_pathway(this_data.C, signal_pathway),
                aliasing_slop=3)
            print("best shift is:",best_shift)
            est_tau.append(best_shift)
            diff = abs(best_shift - programmed_tau)
            diff_tau.append(diff)

