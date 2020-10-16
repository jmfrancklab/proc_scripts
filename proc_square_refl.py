from pyspecdata import *

init_logging("debug")
d = find_file(
    "201009_coilE_1", exp_type="ODNP_NMR_comp/test_equipment", expno="capture1"
)
d.setaxis("ch", r_[1, 2])
d.set_units("t", "s")


class fl_ext(figlist_var):
    def next(self, *arg, **kwargs):
        kwargs.update({"figsize": (9, 6), "legend": True})
        super().next(*arg, **kwargs)

    def abs_re_plot(fl, d, add_text="", show_angle=False):
        for j in range(ndshape(d)["ch"]):
            chlabel = d.getaxis("ch")[j]
            l = fl.plot(
                abs(d["ch", j]),
                linewidth=3,
                alpha=0.5,
                label="CH%d abs " % chlabel + add_text,
            )
            fl.plot(
                d["ch", j].real,
                linewidth=1,
                color=l[0].get_color(),
                alpha=0.5,
                label="CH%d real " % chlabel + add_text,
            )
            fl.plot(
                d["ch", j].imag,
                "--",
                linewidth=1,
                color=l[0].get_color(),
                alpha=0.5,
                label="CH%d imag " % chlabel + add_text,
            )
            if show_angle:
                fl.twinx(orig=False)
                fl.plot(
                    d["ch", j].angle/2/pi,
                    ".",
                    linewidth=1,
                    color=l[0].get_color(),
                    alpha=0.5,
                    label="CH%d angle " % chlabel + add_text,
                )
                fl.twinx(orig=True)
        fl.twinx(orig=False)
        ylabel("phase / cyc", size=10)
        ax = gca()
        gridandtick(ax)
        ax.grid(False)
        fl.twinx(orig=True)
        fl.grid()


with fl_ext() as fl:
    d.ft("t", shift=True)
    d = d["t":(0, 50e6)]
    fl.next("show the frequency distribution")
    forplot = d.C
    forplot[lambda x: abs(x) < 1e-10] = 0
    d["t":(0, 2e6)] = 0
    d["t":(20e6, None)] = 0
    fl.plot(abs(forplot), alpha=0.5, plottype="semilogy", label="orig")
    fl.grid()
    frq_guess = abs(d["ch", 1]).argmax("t").item()
    df = diff(d.getaxis("t")[r_[0, 1]]).item()
    d.ift("t")
    # {{{ determine the frequency from the phase gradient during the pulse
    dt = diff(d.getaxis("t")[r_[0, 1]]).item()
    pulse_slice = d["ch", 0].contiguous(lambda x: abs(x) > 0.5*abs(x).data.max())[0]
    d.setaxis("t", lambda x: x - pulse_slice[0]).register_axis({"t": 0})
    pulse_slice -= pulse_slice[0]
    d = d["t" : tuple(pulse_slice + r_[-0.5e-6, 2e-6])]
    pulse_middle = d["ch", 0]["t" : tuple(pulse_slice + r_[+0.5e-6, -0.5e-6])]
    ph_diff = pulse_middle["t", 1:] / pulse_middle["t", :-1]
    ph_diff.sum("t")
    ph_diff = ph_diff.angle.item()
    frq = ph_diff/dt/2/pi
    # }}}
    print("frq:", frq)
    d *= exp(-1j*2*pi*frq*d.fromaxis("t"))
    ph0 = d["ch", 0].C.sum("t").item()
    ph0 /= abs(ph0)
    d /= ph0
    fl.next("analytic signal -- phase plot", twinx=True)
    for j in range(ndshape(d)["ch"]):
        fl.twinx(orig=True)
        l = fl.plot(abs(d["ch", j]))
        fl.twinx(orig=False)
        fl.plot(d["ch", j].angle, ".", color=l[0].get_color(), alpha=0.1)
    fl.grid()
    fl.next("analytic signal -- abs,re")
    fl.abs_re_plot(d)
    scalar_refl = d["ch", 1]["t":(2e-6, 6e-6)].mean("t").item()
    fl.next("blips")
    blip_range = r_[-0.2e-6, 1.5e-6]
    first_blip = -d["ch", 1:2]["t" : tuple(blip_range)] + scalar_refl
    ph0_blip = first_blip["t", abs(first_blip).argmax("t", raw_index=True).item()]
    ph0_blip /= abs(ph0_blip)
    fl.abs_re_plot(first_blip/ph0_blip, "first")
    secon_blip = d["ch", 1:2]["t" : tuple(blip_range + pulse_slice[1])].setaxis(
        "t", lambda x: x - pulse_slice[1]
    )
    fl.abs_re_plot(secon_blip/ph0_blip, "second", show_angle=True)
