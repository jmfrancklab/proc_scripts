"AG: missing docstring!!!!"
from pylab import *
from pyspecdata import *
from pyspecProcScripts.third_level.analyze_square_refl import analyze_square_refl
from pyspecProcScripts import * 
init_logging("debug")
with fl_mod() as fl:
    for filename, expno, dataset_name in [("210125_sqwv_cap_probe_1", "capture1", "hairpin probe"),
            ('210111_sqwv_sol_probe_1', 'capture1', 'solenoid probe')]:
        logger.info(strm("processing dataset",dataset_name))
        d = find_file(filename,exp_type='ODNP_NMR_comp/test_equip',expno=expno)
        d.set_units('t','s').name('Amplitude').set_units('V')
        d.setaxis("ch", r_[1, 2])
        d.set_units("t", "s")
        analyze_square_refl(d, label=dataset_name, fl=fl,
                show_analytic_signal_phase=False,
                show_analytic_signal_real=False)
