"""
Phasing and Timing Correction Using a Varied Tau Experiment
===========================================================

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
init_logging(level='debug')

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 1
t2, td, vd, power, ph1, ph2 = s.symbols("t2 td vd power ph1 ph2")
f_range = (-400, 400)
filename = '201113_TEMPOL_capillary_probe_var_tau_1'
signal_pathway = {'ph1':1,'ph2':0}
with figlist_var() as fl:
    for nodename,file_location,postproc,label in [
        ('var_tau','ODNP_NMR_comp/test_equipment/var_tau','spincore_var_tau_v1',
            'tau is 1 ms'),
            ]:
        data = find_file(filename,exp_type=file_location,expno=nodename,
                postproc=postproc,lookup=lookup_table)
        data = data['tau',:-7]
        tau_list = list(data.getaxis('tau'))
        data.reorder(['ph1','ph2','tau','t2'])
        data = data['t2':f_range]
        mytable = []
        mytable.append(['programmed tau / ms','estimated tau / ms','difference / ms'])
        for j in range(len(tau_list)):
            tablerow = []
            alias_slop=3
            programmed_tau = tau_list[j]
            tablerow.append(programmed_tau/1e-3)
            logger.info(strm("programmed tau:",programmed_tau))
            this_data = data['tau',j]
            this_data.ift("t2")
            fl.basename = '%0.1f ms'%(programmed_tau/1e-3)
            best_shift = hermitian_function_test(
                select_pathway(this_data, signal_pathway),
                aliasing_slop=alias_slop,
                fl=fl)
            logger.info(strm("best shift is:",best_shift))
            tablerow.append(best_shift/1e-3)
            diff = abs(best_shift - programmed_tau)
            tablerow.append(diff/1e-3)
            mytable.append(tablerow)
        def tabulate(mytable):
            print(' '.join(mytable[0]))
            strlens = [len(j) for j in mytable[0]]
            print(' '.join('-'*j for j in strlens))
            formatstr = ' '.join(f'%{str(j)}.2f' for j in strlens)
            for j in mytable[1:]:
                print(formatstr%tuple(j))
        tabulate(mytable)

