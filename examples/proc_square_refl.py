from pylab import *
from pyspecdata import *
from proc_scripts import *
from proc_scripts.third_level.analyze_square_refl import analyze_square_refl
init_logging("debug")

with fl_mod() as fl:
    for filename, expno, dataset_name in [("201218_sqwv_cap_probe_2", "capture1", "hairpin probe"),
            ('201228_sqwv_sol_probe_1', 'capture1', 'solenoid probe')]:
        print("processing dataset",dataset_name)
        d = find_file(filename,exp_type='ODNP_NMR_comp/test_equip',expno=expno)
        d.set_units('t','s').name('Amplitude').set_units('V')
        d.setaxis("ch", r_[1, 2])
        d.set_units("t", "s")
        analyze_square_refl(d, label=dataset_name, fl=fl,
                show_analytic_signal_phase=False,
                show_analytic_signal_real=False)
