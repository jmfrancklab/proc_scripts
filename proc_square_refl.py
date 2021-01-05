from pyspecdata import *
from scipy.signal import tukey
from scipy.optimize import minimize,leastsq
from sympy import symbols
init_logging("debug")

d = find_file('201218_sqwv_cap_probe_1',exp_type='ODNP_NMR_comp/test_equip',expno='capture1')
s = find_file(
    "201228_sqwv_sol_probe_1", exp_type="ODNP_NMR_comp/test_equip", expno="capture1"
)
d.set_units('t','s').name('Amplitude').set_units('V')
s.set_units('t','s').name('Amplitude').set_units('V')
d.setaxis("ch", r_[1, 2])
s.setaxis('ch',r_[1,2])
d.set_units("t", "s")
s.set_units("t","s")
with figlist_var() as fl:
    fl.next('Raw signal')
    fl.plot(s['ch',0], alpha=0.5, label='solenoid control')
    fl.plot(s['ch',1], alpha=0.5, label='solenoid reflection')
    fl.plot(d['ch',0], alpha=0.5, label='capillary control')
    fl.plot(d['ch',1], alpha=0.5, label='capillary reflection')
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
    d.ft("t")#, shift=True)
    s.ft('t')
    d = 2*d['t':(0,None)]
    s = 2*s['t':(0,None)]
    #in the previous version we had multiplied data by 2 because the euation 1/2a*exp(iwt)+aexp(-iwt) and the 2
    #negated the half.... this is not done here. before we did it because the demodulation looked weird w/o it.
    #I am assuming this is not needed here as the demodulation either won't be used or is not affected anymore.
    d = d["t":(0, 50e6)]
    s = s['t':(0,50e6)]
    fl.next("show the frequency distribution")
    forplot = d.C
    forplot_s = s.C
    forplot[lambda x: abs(x) < 1e-10] = 0 #filter low frequency noise out
    forplot_s[lambda x:abs(x) < 1e-10] = 0
    frq_guess = abs(d["ch", 1]).argmax("t").item() #this is to find the peak? argmax on reflection ch would be the peaks
    frq_guess_s = abs(s['ch',1]).argmax('t').item()
    frq_range = r_[-30,30]*1e6 + frq_guess
    frq_range_s = r_[-30,30]*1e6 + frq_guess_s
    d["t":(0, frq_range[0])] = 0 #again filtering out noise outside of blips
    s['t':(0,frq_range_s[0])] = 0
    d["t":(frq_range[1], None)] = 0
    s['t':(frq_range_s[1],None)] = 0
    # {{{ shouldn't have to do it this way, but something weird going on w/ aligndata
    tukey_filter = d.fromaxis("t")["t":tuple(frq_range)].run(lambda x: tukey(len(x)))
    tukey_filter_s = s.fromaxis('t')['t':tuple(frq_range_s)].run(lambda x: tukey(len(x)))
    d["t":tuple(frq_range)] *= tukey_filter
    s['t':tuple(frq_range_s)] *= tukey_filter_s
    for j in d.getaxis('ch'):
        fl.plot(abs(forplot)['ch':j], alpha=0.5, plottype="semilogy", label=f"CH{j} orig capillary")
        fl.plot(abs(d)['ch':j][lambda x: abs(x) > 1e-10], alpha=0.5, plottype="semilogy", label=f"CH{j} filtered capillary") 
    for j in s.getaxis('ch'):
        fl.plot(abs(forplot_s)['ch':j], alpha=0.5, plottype='semilogy', label=f'CH{j} original solenoid')
        fl.plot(abs(s)['ch':j][lambda x: abs(x) > 1e-10], alpha=0.5, plottype='semilogy', label=f'CH{j} filtered solenoid')
    fl.grid()
    df = diff(d.getaxis("t")[r_[0, 1]]).item()
    sf = diff(s.getaxis('t')[r_[0, 1]]).item()
    d.ift("t")
    s.ift('t')
    # {{{ determine the frequency from the phase gradient during the pulse
    dt = diff(d.getaxis("t")[r_[0, 1]]).item()
    st = diff(s.getaxis('t')[r_[0, 1]]).item()
    pulse_slice = d["ch", 0].contiguous(lambda x: abs(x) > 0.5*abs(x).data.max())[0] #defines pulse slice based on control signal
    pulse_slice_s = s['ch',0].contiguous(lambda x: abs(x) > 0.5*abs(x).data.max())[0]
    d.setaxis("t", lambda x: x - pulse_slice[0]).register_axis({"t": 0}) #resets t axis around pulse slice
    s.setaxis('t', lambda x: x - pulse_slice_s[0]).register_axis({'t':0})
    pulse_slice -= pulse_slice[0]
    pulse_slice_s -= pulse_slice_s[0]
    #{{{ Not sure what this portion is doing...is this similar to a time shift? or first order phase correction?
    #fl.show();quit()
    d = d["t" : tuple(pulse_slice + r_[-0.5e-6, 5e-6])]
    s = s['t':tuple(pulse_slice_s + r_[-0.5e-6, 5e-6])]
    pulse_middle = d["ch", 0]["t" : tuple(pulse_slice + r_[+0.5e-6, -0.5e-6])]
    pulse_middle_s = s['ch',0]['t':tuple(pulse_slice_s + r_[+0.5e-6,-0.5e-6])]
    ph_diff = pulse_middle["t", 1:] / pulse_middle["t", :-1]
    ph_diff_s = pulse_middle_s['t',1:] / pulse_middle_s['t', :-1]
    ph_diff.sum("t")
    ph_diff_s.sum('t')
    ph_diff = ph_diff.angle.item()
    ph_diff_s = ph_diff_s.angle.item()
    frq = ph_diff/dt/2/pi
    frq_s = ph_diff_s/st/2/pi
    #}}}
    # }}}
    print("frq:", frq)
    d *= exp(-1j*2*pi*frq*d.fromaxis("t")) #convolution
    s *= exp(-1j*2*pi*frq_s*s.fromaxis('t'))
    ph0 = d["ch", 0].C.sum("t").item() #pseudo 0th order phase correction
    ph0_s = s['ch',0].C.sum('t').item()
    ph0 /= abs(ph0)
    ph0_s /= abs(ph0_s)
    d /= ph0
    s /= ph0_s
    fl.next("analytic signal -- phase plot", twinx=True)
    for j in range(ndshape(d)["ch"]):
        fl.twinx(orig=True)
        l = fl.plot(abs(d["ch", j]))
        fl.twinx(orig=False)
        fl.plot(d["ch", j].angle, ".", color=l[0].get_color(), alpha=0.1,label='capillary probe angle')
    for j in range(ndshape(s)['ch']):
        fl.twinx(orig=True)
        ls = fl.plot(abs(s['ch',j]))
        fl.twinx(orig=False)
        fl.plot(s['ch',j].angle,".",color=ls[0].get_color(),alpha=0.1,label='solenoid probe angle')
    fl.grid()
    fl.next("analytic signal -- abs,re")
    fl.abs_re_plot(d,add_text="capillary probe")
    fl.abs_re_plot(s,add_text="solenoid probe")
    scalar_refl = d["ch", 1]["t":(1e-6, 2e-6)].mean("t").item()
    scalar_refl = s['ch',1]['t':(1e-6,1.2e-6)].mean('t').item()
    fl.show();quit()
    fl.next("blips")
    blip_range = r_[-0.1e-6, 2e-6] #defining decay slice
    first_blip = -d["ch", 1:2]["t" : tuple(blip_range)] + scalar_refl #correcting first blip
    #{{{ doing 0th order correction type thing again? why? we did this in lines 99-101...
    ph0_blip = first_blip["t", abs(first_blip).argmax("t", raw_index=True).item()]
    ph0_blip /= abs(ph0_blip)
    fl.abs_re_plot(first_blip/ph0_blip, "first")
    secon_blip = d["ch", 1:2]["t" : tuple(blip_range + pulse_slice[1])].setaxis(
        "t", lambda x: x - pulse_slice[1]
    )
    fl.abs_re_plot(secon_blip/ph0_blip, "second", show_angle=True)
    decay = (abs(first_blip)+abs(secon_blip))/2
    fl.next('decay')
    fl.plot(decay)
    #fl.show();quit()
    #decay_start = decay.argmax('t').item()
    #decay = decay['t':(decay_start,None)]
    decay = decay['t':(41e-9,1200)]
    fl.next('Plotting the decay slice')
    fl.plot(decay, linewidth=3, alpha=0.3, color='k')
    print(decay.getaxis('ch'))
    decay = decay['ch',0]
    print(ndshape(decay))
    f = fitdata(decay)
    A,B,C,t = symbols("A B C t",real=True)
    f.functional_form = A*e**(-t*B)
    fl.next('fit for capillary probe')
    fl.plot(decay,'o',label='data')
    f.fit()
    f.set_units('t','ns')
    print("output:",f.output())
    print("latex:",f.latex())
    Q = 1./f.output('B')*2*pi*14893772
    fl.plot(f.eval(100).set_units('t','s'),label='fit, Q=%0.1f'%Q)
quit()

