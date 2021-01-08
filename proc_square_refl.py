from pylab import *
from pyspecdata import *
from scipy.signal import tukey
from scipy.optimize import minimize,leastsq
import sympy as s
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
                fl.twinx(orig=True)
        fl.twinx(orig=False)
        ylabel("phase / cyc", size=10)
        ax = gca()
        gridandtick(ax)
        ax.grid(False)
        fl.twinx(orig=True)
        fl.grid()

filename, expno = ["201228_sqwv_sol_probe_1", "capture1"]
#filename, expno = ['201218_sqwv_cap_probe_1', 'capture1']
frq_bw = 15e6
keep_after_pulse = 2e-6 # amount of the time axis after the pulse we
#                         want hanging around -- this can change
#                         depending on the Q
blip_range = r_[-0.1e-6, keep_after_pulse] # where the blip is relative to the edge of the pulse
dataset_name = 'capillary'
d = find_file(filename,exp_type='ODNP_NMR_comp/test_equip',expno=expno)
d.set_units('t','s').name('Amplitude').set_units('V')
d.setaxis("ch", r_[1, 2])
d.set_units("t", "s")


with fl_ext() as fl:
    d.ft("t", shift=True)
    d = d['t':(-50e6,50e6)] # slice out a reasonable range
    d = 2*d['t':(0,None)] # multiply data by 2 because the equation
    #                       1/2a*exp(iwt)+aexp(-iwt) and the 2 negated the
    #                       half.
    fl.next("show the frequency distribution")
    for j in d.getaxis('ch'):
        fl.plot(abs(d)['ch':j][lambda x: abs(x) > 1e-10], alpha=0.5,
                plottype="semilogy", label=f"CH{j} {dataset_name}") 
    frq_guess = abs(d["ch", 0]).argmax("t").item() # find the peak
    frq_range = r_[-frq_bw,frq_bw]*0.5 + frq_guess
    d["t":(0, frq_range[0])] = 0 # throw out everything after the top of our slice
    d["t":(frq_range[1], None)] = 0 # set everything else to zero
    # {{{ shouldn't have to do it this way, but something weird going on w/ aligndata
    tukey_filter = d.fromaxis("t")["t":tuple(frq_range)].run(lambda x: tukey(len(x)))
    d["t":tuple(frq_range)] *= tukey_filter
    for j in d.getaxis('ch'):
        fl.plot(abs(d)['ch':j][lambda x: abs(x) > 1e-10], alpha=0.5,
                plottype="semilogy", label=f"CH{j} filtered {dataset_name}") 
    fl.grid()
    df = diff(d.getaxis("t")[r_[0, 1]]).item()
    d.ift("t")
    # {{{ slice out pulse and surrounding and start time axis with pulse
    dt = diff(d.getaxis("t")[r_[0, 1]]).item()
    pulse_slice = d["ch", 0].contiguous(lambda x:
            abs(x) > 0.5*abs(x).data.max())[0] # defines pulse slice based on control signal
    d.setaxis("t", lambda x: x - pulse_slice[0]).register_axis({"t": 0}) # set t=0 to the beginning of the pulse
    pulse_slice -= pulse_slice[0]
    d = d["t" : tuple(pulse_slice + r_[-0.5e-6, keep_after_pulse])] # take a slice that starts 0.5 μs before the pulse and ends 2 μs after the pulse
    # }}}
    # {{{ determine ∂φ/∂t during the pulse, use to det frq
    pulse_middle_slice = tuple(pulse_slice + r_[+0.5e-6, -0.5e-6])
    pulse_middle = d["ch", 0]["t": pulse_middle_slice]
    ph_diff = pulse_middle["t", 1:] / pulse_middle["t", :-1]
    ph_diff *= abs(pulse_middle["t", :-1])
    # at this point, we have something of the form
    # ρe^{iΔφ} for each point -- average them
    ph_diff.sum("t")
    ph_diff = ph_diff.angle.item() # Δφ above
    frq = ph_diff/dt/2/pi
    # }}}
    print("frq:", frq)
    d *= exp(-1j*2*pi*frq*d.fromaxis("t")) # mix down
    d.ft('t')
    fl.next('after slice and mix down freq domain')
    fl.plot(abs(d), label=dataset_name)
    d.ift('t')
    # {{{ zeroth order phase correction
    ph0 = d["ch", 0]['t':pulse_middle_slice].C.sum("t").item()
    ph0 /= abs(ph0)
    d /= ph0
    # }}}
    fl.next("analytic signal", twinx=True)
    fl.complex_plot(d, label=dataset_name, show_phase=True)
    fl.grid()
    scalar_refl = d["ch", 1]["t":(keep_after_pulse, pulse_middle_slice[-1])].mean("t").item()
    fl.next("blips")
    first_blip = -d["ch", 1:2]["t":tuple(blip_range)] + scalar_refl # correcting first blip
    fl.complex_plot(first_blip, "first", show_phase=True, show_real=False)
    secon_blip = d["ch", 1:2]["t" : tuple(blip_range + pulse_slice[1])].setaxis(
        "t", lambda x: x - pulse_slice[1]
    )
    fl.complex_plot(secon_blip, "second", show_phase=True, show_real=False)
    decay = abs(secon_blip)['ch',0] # we need the ch axis for the complex plot,
    #                                 but it complicates things now
    decay_start = decay.C.argmax('t').item()
    decay = decay['t':(decay_start,None)]
    f = fitdata(decay)
    A,B,C,t = s.symbols("A B C t",real=True)
    f.functional_form = A*s.exp(-t*B)+C
    f.fit()
    print("output:",f.output())
    print("latex:",f.latex())
    Q = 1./f.output('B')*2*pi*frq
    fl.plot(f.eval(100).set_units('t','s'),'k--', alpha=0.8, label='fit, Q=%0.1f'%Q)
