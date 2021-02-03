import pylab as plb
from pylab import r_, pi
import pyspecdata as psp
from scipy.signal import tukey
from scipy.optimize import minimize,leastsq
import sympy as s
import matplotlib.pyplot as plt
def analyze_square_refl(d, label='', fl=None,
        frq_bw=15e6,
        keep_after_pulse=2e-6, # amount of the time
        #                        axis after the pulse we want hanging around --
        #                        this can change depending on the Q
        blip_range = [-0.1e-6,None],
        show_analytic_signal_phase = True,
        show_analytic_signal_real = False,
        ):
    r"""
    Parameters
    ==========
    fl: figlist_var child class
        In addition to standard figlist_var methods, this must also include a
        `complex_plot` method that returns list "colors"
    """
    if len(label)>0:
        fl.basename = label
    orig_basename = fl.basename
    if blip_range[-1] is None:
        blip_range[-1] = keep_after_pulse
    d.ft("t", shift=True)
    d.ift('t')
    d.ft('t')
    d = d['t':(-50e6,50e6)] # slice out a reasonable range
    d['t':(None,0)]['t',:-1] = 0
    d *= 2 # multiply data by 2 because the equation
    #                       1/2a*exp(iwt)+aexp(-iwt) and the 2 negated the
    #                       half.
    if fl is not None:
        fl.next("show the frequency distribution")
        for j in d.getaxis('ch'):
            fl.plot(abs(d)['ch':j][lambda x: abs(x) > 1e-10], alpha=0.5,
                    plottype="semilogy", label=f"CH{j} {label}") 
    frq_guess = abs(d["ch", 0]).argmax("t").item() # find the peak
    frq_range = r_[-frq_bw, frq_bw] * 0.5 + frq_guess
    d["t":(0, frq_range[0])] = 0 # set everything else to zero
    d["t":(frq_range[1], None)] = 0
    # {{{ shouldn't have to do it this way, but something weird going on w/ aligndata
    tukey_filter = d.fromaxis("t").C
    tukey_filter = tukey_filter["t":tuple(frq_range)].run(lambda x: tukey(len(x)))
    d["t":tuple(frq_range)] *= tukey_filter
    if fl is not None:
        for j in d.getaxis('ch'):
            fl.plot(abs(d)['ch':j][lambda x: abs(x) > 1e-10], alpha=0.5,
                    plottype="semilogy", label=f"CH{j} filtered {label}") 
        fl.grid()
    df = plb.diff(d.getaxis("t")[r_[0, 1]]).item()
    d.ift("t")
    # {{{ slice out pulse and surrounding and start time axis with pulse
    dt = plb.diff(d.getaxis("t")[r_[0, 1]]).item()
    pulse_slice = d["ch", 0].contiguous(lambda x:
                                        abs(x) > 0.5 * abs(x).data.max())[0] # defines pulse slice based on control signal
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
    frq = ph_diff / dt / 2 / pi
    # }}}
    print("frq:", frq)
    d *= plb.exp(-1j * 2 * pi * frq * d.fromaxis("t")) # mix down
    d.ft('t')
    if fl is not None: fl.next('after slice and mix down freq domain')
    if fl is not None: fl.plot(abs(d), label=label)
    d.ift('t')
    # {{{ zeroth order phase correction
    ph0 = d["ch", 0]['t':pulse_middle_slice].C.mean("t").item()
    pulse_middle_amp = abs(ph0)
    ph0 /= pulse_middle_amp
    d /= ph0
    # }}}
    if fl is not None:
        fl.basename = None
        fl.next("analytic signal")
        colors = fl.complex_plot(d, label=label,
                show_phase=show_analytic_signal_phase,
                show_real=show_analytic_signal_real,alpha=0.5)
        # {{{ print the carrier
        # transform goes to "display", which is pixels
        # "inverted" goes back
        # here, I could set a mixed transformation,
        # but I think it's more understandable to just do manually
        ax = plb.gca()
        print("the amplitude is",pulse_middle_amp)
        _,y = ax.transData.transform(r_[0.0, pulse_middle_amp])
        fontsize = 16
        nfigures = len(fl.figurelist)
        y -= fontsize*(nfigures/3)
        _,y = ax.transAxes.inverted().transform(r_[0, y])
        plb.text(
                x=0.5,
                y=y,
                s=r'$\nu_{Tx}=%0.6f$ MHz'%(frq/1e6),
                va='top',
                ha='center',
                size=fontsize,
                transform=ax.transAxes, # display
                color=colors[0],
                )
        # }}}
        fl.basename = orig_basename
    scalar_refl = d["ch", 1]["t":(keep_after_pulse, pulse_middle_slice[-1])].mean("t").item()
    if label == "hairpin probe":
        color='darkorange'
    if label == 'solenoid probe':
        color='red'
    if fl is not None:
        fl.basename=None
        fl.next("blips")
        first_blip = -d["ch", 1:2]["t":tuple(blip_range)] + scalar_refl # correcting first blip
        fl.complex_plot(first_blip, "first", show_phase=False, show_real=False,alpha=0.2,linestyle="--",linewidth=1,color=color)
        secon_blip = d["ch", 1:2]["t" : tuple(blip_range + pulse_slice[1])].setaxis(
            "t", lambda x: x - pulse_slice[1]
        )
        fl.complex_plot(secon_blip, "second", show_phase=True, show_real=False,alpha=0.8,color=color)
    secon_blip = secon_blip['ch', 0] # we need the ch axis for the complex plot,
    #                                  but it complicates things now
    decay = abs(secon_blip)
    decay_start = decay.C.argmax('t').item()
    decay = decay['t':(decay_start,None)]
    f = psp.fitdata(decay)
    A,R,C,t = s.symbols("A R C t",real=True)
    f.functional_form = A*s.exp(-t*R)+C
    f.set_guess({A:0.3, R: 1 / 20 * 2 * pi * frq})
    f.fit()
    print("output:",f.output())
    print("latex:",f.latex())
    # {{{ calculate frequency offset
    decay_timescale = 3./f.output('R')
    dt = plb.diff(decay.getaxis("t")[r_[0, 1]]).item()
    phases = secon_blip['t':(150*10**-9,((54. / f.output('R') * 2 * pi * frq)*10**-9))]
    print(phases)
    quit()
    frq_offset = (phases['t',1:]/phases['t',:-1]*abs(phases['t',:-1])
            ).sum('t').angle.item()/dt/2/pi
    if fl is not None:
        ax = fl.twinx(orig=False)
        x = plb.mean(phases.getaxis('t'))/1e-9
        y = plb.mean(phases.angle.data)/1.25/pi
        ax.text(
                x=x,
                y=y,
                s=' '*5+r'$\Delta\nu=%0.3g$ kHz'%(frq_offset/1e3),
                va='bottom',
                ha='left',
                size=fontsize,
                transform=ax.transData,
                color=color,
                )
        fl.twinx(orig=True)
    print('frq_offset',frq_offset,"for",label)
    # }}}
    Q = 1. / f.output('R') * 2 * pi * frq
    if fl is not None:
        ax = fl.twinx(orig=False)
        x = plb.mean(phases.getaxis('t'))/1e-9
        y = plb.mean(phases.angle.data)/6/pi
        ax.text(
                x=x,
                y=y,
                s=' '*5+r'$Q=%d$'%(Q),
                va='bottom',
                ha='left',
                size=fontsize,
                transform=ax.transData,
                color=color,
                )
        fl.twinx(orig=True)

    if fl is not None: fl.plot(f.eval(100).set_units('t','s'),'k--', alpha=0.8)#, label='fit, Q=%0.1f'%Q)

