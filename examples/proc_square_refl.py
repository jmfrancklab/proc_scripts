from pylab import *
from pyspecdata import *
from proc_scripts.third_level.analyze_square_refl import analyze_square_refl
init_logging("debug")

class fl_ext(figlist_var):
    def next(self, *arg, **kwargs):
        kwargs.update({"figsize": (9, 6), "legend": True})
        super().next(*arg, **kwargs)

    def complex_plot(fl, d, label="", show_phase=False, show_real=True):
        colors = []
        for j in range(ndshape(d)["ch"]):
            chlabel = d.getaxis("ch")[j]
            l = fl.plot(
                abs(d["ch", j]),
                linewidth=3,
                alpha=0.5,
                label="CH%d abs " % chlabel + label,
            )
            colors.append(l[0].get_color())
            if show_real:
                fl.plot(
                    d["ch", j].real,
                    linewidth=1,
                    color=colors[-1],
                    alpha=0.5,
                    label="CH%d real " % chlabel + label,
                )
                fl.plot(
                    d["ch", j].imag,
                    "--",
                    linewidth=1,
                    color=colors[-1],
                    alpha=0.5,
                    label="CH%d imag " % chlabel + label,
                )
            fl.grid()
            if show_phase:
                fl.twinx(orig=False)
                fl.plot(
                    d["ch", j].angle/2/pi,
                    ".",
                    linewidth=1,
                    color=colors[-1],
                    alpha=0.3,
                    label="CH%d angle " % chlabel + label,
                )
                ylabel("phase / cyc", size=10)
                ax2=gca()
                gridandtick(ax2, use_grid=False)
                ax2.grid(False)
                fl.twinx(orig=True)
        return colors

with fl_ext() as fl:
    for filename, expno, dataset_name in [("210125_sqwv_cap_probe_1", "capture1", "hairpin probe"),
            ('210111_sqwv_sol_probe_1', 'capture1', 'solenoid probe')]:
        print("processing dataset",dataset_name)
        d = find_file(filename,exp_type='ODNP_NMR_comp/test_equip',expno=expno)
        d.set_units('t','s').name('Amplitude').set_units('V')
        d.setaxis("ch", r_[1, 2])
        d.set_units("t", "s")
        analyze_square_refl(d, label=dataset_name, fl=fl,
                show_analytic_signal_phase=False,
                show_analytic_signal_real=False)
