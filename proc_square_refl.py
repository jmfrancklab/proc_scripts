from pylab import *
from pyspecdata import *
from proc_scripts.third_level.analyze_square_refl import analyze_square_refl
init_logging("debug")

class fl_ext(figlist_var):
    def next(self, *arg, **kwargs):
        kwargs.update({"figsize": (9, 6), "legend": True})
        super().next(*arg, **kwargs)

    def complex_plot(fl, d, label="", show_phase=False, show_real=True):
        for j in range(ndshape(d)["ch"]):
            if show_phase:
                fl.twinx(orig=True)
            chlabel = d.getaxis("ch")[j]
            l = fl.plot(
                abs(d["ch", j]),
                linewidth=3,
                alpha=0.5,
                label="CH%d abs " % chlabel + label,
            )
            if show_real:
                fl.plot(
                    d["ch", j].real,
                    linewidth=1,
                    color=l[0].get_color(),
                    alpha=0.5,
                    label="CH%d real " % chlabel + label,
                )
                fl.plot(
                    d["ch", j].imag,
                    "--",
                    linewidth=1,
                    color=l[0].get_color(),
                    alpha=0.5,
                    label="CH%d imag " % chlabel + label,
                )
            if show_phase:
                fl.twinx(orig=False)
                fl.plot(
                    d["ch", j].angle/2/pi,
                    ".",
                    linewidth=1,
                    color=l[0].get_color(),
                    alpha=0.3,
                    label="CH%d angle " % chlabel + label,
                )
                ylabel("phase / cyc", size=10)
                fl.twinx(orig=True)
            fl.grid()

filename, expno, dataset_name = ["201228_sqwv_sol_probe_1", "capture1", "capillary"]
#filename, expno = ['201218_sqwv_cap_probe_1', 'capture1']
d = find_file(filename,exp_type='ODNP_NMR_comp/test_equip',expno=expno)
d.set_units('t','s').name('Amplitude').set_units('V')
d.setaxis("ch", r_[1, 2])
d.set_units("t", "s")

with fl_ext() as fl:
    analyze_square_refl(d, label=dataset_name, fl=fl)
