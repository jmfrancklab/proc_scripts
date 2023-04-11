"""
Quantitative ESR
================
In this example, we baseline correct and
integrate ESR data, and compare the result to
the result of integration inside XEPR.

This makes use of the package `pint`, which is a
very nice package for handling units.
"""
from pyspecdata import *
import numpy as np
from matplotlib.pyplot import axvline, axhline, gca
from pint import UnitRegistry
from pyspecdata.datadir import pyspec_config # piggyback on _pyspecdata
def act_conc(filename,
        background,
        thislabel,
        exp_type='francklab_esr/alex',
        pushout=3,
        fl=None):
    default_Q = float(pyspec_config.get_setting("default Q",default=4700))
    dint_conversion = 1.66e2
    if dint_conversion is None:
        raise ValueError("You cannot run this file without first running a file that calibrates the double integral to concentration conversion!")
    dint_conversion = float(dint_conversion)
    ureg = UnitRegistry(system="mks",autoconvert_offset_to_baseunit=True, auto_reduce_dimensions=True)
    Q_ = ureg.Quantity
    colors = plt.rcParams["axes.prop_cycle"]() # this is the default matplotlib cycler for line styles
    fieldaxis = '$B_0$'
    myconcs = []
    with figlist_var() as fl:
        background = find_file(background,exp_type = "francklab_esr/alex")['harmonic',0]
        background -= background[fieldaxis, -100:].data.mean()
        d = find_file(filename, exp_type="francklab_esr/alex")['harmonic',0]
        G_R = Q_(*d.get_prop("Gain"))
        C_t = Q_(*d.get_prop("ConvTime"))
        power =  Q_(*d.get_prop("Power"))
        B_m = Q_(*d.get_prop("ModAmp"))
        Q = Q_(default_Q,'dimensionless') # hard set Q value
        n_B = Q_(1,'dimensionless') # calculate this
        S = Q_(0.5,'dimensionless')
        c = Q_(1,'dimensionless') # the first fraction on pg 2-17 -- essentially the conversion factor
        signal_denom = G_R * C_t * sqrt(power) * B_m * n_B * S * (S+1) * Q
        signal_denom = signal_denom.to(Q_('G')*sqrt(Q_('W'))*Q_('s'))
        d.set_plot_color(next(colors)['color'])
        forplot = background.C.integrate(fieldaxis, cumulative=True)
        if fl is not None:
            fl.next('show background subtraction in abs mode')
            forplot.set_plot_color(d.get_plot_color())
            fl.plot(forplot, "k", alpha=0.2)
        d -= d[fieldaxis, -100:].data.mean()
        forplot = d.C.integrate(fieldaxis, cumulative = True)
        if fl is not None:
            fl.plot(forplot,alpha=0.5, label=thislabel)
        d -= background
        d_abs = d.C.integrate(fieldaxis, cumulative=True)
        if fl is not None:
            fl.next("absorption, bg. no bl.")
            fl.plot(d_abs, alpha=0.5)
        peaklist = d_abs.contiguous(lambda x: abs(x) > abs(x).data.max()/2)[:3,:]
        peaksize = np.diff(peaklist, axis=1).mean()
        specrange = (peaklist.ravel().min(),peaklist.ravel().max())
        d_int_direct = d_abs.C
        if fl is not None:
            fl.next("d_abs. int, bs. no bl.")
            fl.plot(d_int_direct.integrate(fieldaxis, cumulative=True), alpha=0.5)
        generous_limits = specrange+r_[-pushout*peaksize,+pushout*peaksize]
        if fl is not None:
            fl.next("absorption, bg. no bl.")
            for j in r_[np.array(specrange)]:
                axvline(x=j/1e3, alpha=0.1, color='k', ls=':')
            for j in r_[generous_limits]:
                axvline(x=j/1e3, alpha=0.1, color='k', ls='-')
        d_baseline = d_abs[fieldaxis, lambda x: np.logical_or(x<generous_limits[0], x>generous_limits[1])]
        if fl is not None:
            fl.next("for baseline")
            fl.plot(d_baseline, '.', alpha=0.3, human_units=False)
        middle_field = np.diff(d_baseline.getaxis(fieldaxis)[r_[0,-1]]).item()
        d_baseline.setaxis(fieldaxis, lambda x: x-middle_field)
        c = d_baseline.polyfit(fieldaxis, order=7)
        polybaseline = d.fromaxis(fieldaxis).setaxis(fieldaxis, lambda x: x-middle_field).eval_poly(c,fieldaxis)
        polybaseline.setaxis(fieldaxis, lambda x: x+middle_field)
        if fl is not None:
            fl.plot(polybaseline, alpha=0.5, human_units=False)
            fl.next('absorption, bg. no bl.')
            fl.plot(polybaseline, alpha = 0.5)
        d_abs -= polybaseline
        if fl is not None:
            fl.next('absorption, with baseline')
            fl.plot(d_abs, alpha=0.5)
        d_abs.integrate(fieldaxis, cumulative=True)
        d_abs /= signal_denom.magnitude
        d_abs *= dint_conversion
        d_abs /= 1e-3
        aver = d_abs[fieldaxis:(3.57e3,None)].mean(fieldaxis).item()
        print("Average concentration:",aver)
        d_abs.name('conc').set_units('mM')
        if fl is not None:
            fl.next("dblint ÷ denom * conversion")
            fl.plot(d_abs, alpha=0.5, label=f"{thislabel}, %0.4f mM"%aver)
            fl.grid()
        expected_QESR = (d_abs).data[-100:].mean()
        return expected_QESR 
